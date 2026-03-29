#include "sim_manager.hh"
#include <cstring>
#include <algorithm>
#include <random>

sim_manager::~sim_manager() {
	for (int i = 0; i < n_obj; i++) delete objs[i];
    if (avg_velx != NULL) delete[] avg_velx;
    if (avg_vely != NULL) delete[] avg_vely;
    if (avg_velz != NULL) delete[] avg_velz;
    if (avg_velN != NULL) delete[] avg_velN;
    if (avg_x != NULL) delete[] avg_x;
    if (avg_y != NULL) delete[] avg_y;
    if (avg_z != NULL) delete[] avg_z;
    if (avg_N != NULL) delete[] avg_N;
	if (objs != NULL) delete[] objs;
    if(sm_array!=NULL) delete [] sm_array;
}

double sim_manager::cfl(double dx) {

	// the minimum time scale allowed by shear wave of the solid
	double xxsp3 = 3.0/(dx*dx);
	double dt_s = 100000;
	for (int i = 0; i < n_obj; i++) {
		double odt = sm_array[i].cfl(dx);
		if (odt < dt_s) dt_s = odt;
	}
	// the time scale set by fluid viscosity
	double dt_f = 0.5* fm.rho/(fm.mu*xxsp3)*fm.dt_pad;

	return (dt_s<dt_f)?dt_s:dt_f;
}

/** Set up
 * CFL max dt and transition zone width.
 */
void sim_manager::setup_const(){
	double dh=spars->min_dh();
	// max dt allowed by cfl condition
	dt_reg = cfl(dh);
	// transition zone widht
	eps = wt_n*dh;
	eps_inv = 1./eps;
	tf1 =0.5/eps;
	tf2 =0.5/M_PI;
	tf3 =M_PI/eps;

    // User provided time step could be smaller
    if(spars->dt < dt_reg) dt_reg = spars->dt;

    // Set the anchoring acceleration
    K_stiff = 0.0025/(dt_reg*dt_reg);
}

/** This is to check for advective CFL condition after fluid velocity has been initialized */
void sim_manager::obligatory_cfl_recheck(){
        double dh=spars->min_dh();
        double dt_min = fm.dt_pad*dh/vmax;
        if(dt_reg>dt_min) {
            dt_reg = dt_min;
            // User provided time step could be smaller
            if(spars->dt < dt_reg) dt_reg = spars->dt;
        }
        // Set the anchoring acceleration
        K_stiff = 0.0025/(dt_reg*dt_reg);
}

void sim_manager::create_objs(){
    for(int i=0;i<n_obj;i++){
        sm_array[i] = sl_mat(spars, i);
        // HACK HACK HACK:
        sm_array[i].ex_mu = sm_array[i].ex_visc_mult*fm.mu;

        double *basics = spars->basic_specs + n_basic_specs*i;
        double *extras = spars->extra_specs + n_extra_specs*i;

        int obj_type_index = spars->object_list[i];

        objs [i] = object::alloc(obj_type_index, basics, extras);
    }
}

void sim_manager::setup_isotropic_turbulence() {
    if(!turb_iso_enable) return;
    turb_modes.clear();
    turb_modes.reserve(turb_mode_count);
    std::mt19937 gen(static_cast<unsigned int>(turb_seed));
    std::uniform_int_distribution<int> kdist(-turb_kf2, turb_kf2);
    std::uniform_real_distribution<double> pdist(0.0, 2.0*M_PI);
    std::uniform_real_distribution<double> udist(-1.0, 1.0);

    const double nu = fm.mu / fm.rho;
    turb_target_energy = turb_re_lam*std::pow(nu*turb_kd, 2)/std::sqrt(20./3.);

    int guard = 0;
    while(static_cast<int>(turb_modes.size()) < turb_mode_count && guard < 200*turb_mode_count) {
        guard++;
        int kx = kdist(gen), ky = kdist(gen), kz = kdist(gen);
        if(kx==0 && ky==0 && kz==0) continue;
        int k2 = kx*kx + ky*ky + kz*kz;
        if(k2 > turb_kf2*turb_kf2) continue;
        double kmag = std::sqrt(static_cast<double>(k2));
        double kvec[3] = {static_cast<double>(kx), static_cast<double>(ky), static_cast<double>(kz)};

        // build two random orthonormal vectors in plane orthogonal to k
        double r0 = udist(gen), r1 = udist(gen), r2 = udist(gen);
        double dot = r0*kvec[0] + r1*kvec[1] + r2*kvec[2];
        r0 -= dot*kvec[0]/k2;
        r1 -= dot*kvec[1]/k2;
        r2 -= dot*kvec[2]/k2;
        double nr = std::sqrt(r0*r0 + r1*r1 + r2*r2);
        if(nr < 1e-12) continue;
        r0 /= nr; r1 /= nr; r2 /= nr;
        double s0 = (kvec[1]*r2 - kvec[2]*r1)/kmag;
        double s1 = (kvec[2]*r0 - kvec[0]*r2)/kmag;
        double s2 = (kvec[0]*r1 - kvec[1]*r0)/kmag;

        double phase1 = pdist(gen), phase2 = pdist(gen);
        // Rogallo-like shape: low-k ~ k, high-k ~ k^(-5/3)
        double ek0 = std::sqrt(9./11./std::max(1, turb_kf2) * static_cast<double>(k2)/(turb_kf2*turb_kf2));
        double ek1 = std::sqrt(9./11./std::max(1, turb_kf2) * std::pow(std::max(kmag, 1e-12)/std::max(1, turb_kf2), -5./3.));
        double ek = (k2 <= turb_kf2*turb_kf2) ? ek0 : ek1;
        double amp = ek/std::sqrt(std::max(1, turb_mode_count));
        turb_mode mode;
        mode.k[0]=kx; mode.k[1]=ky; mode.k[2]=kz;
        mode.a[0]=amp*(std::cos(phase1)*r0 + std::sin(phase1)*s0);
        mode.a[1]=amp*(std::cos(phase1)*r1 + std::sin(phase1)*s1);
        mode.a[2]=amp*(std::cos(phase1)*r2 + std::sin(phase1)*s2);
        mode.b[0]=amp*(std::cos(phase2)*r0 + std::sin(phase2)*s0);
        mode.b[1]=amp*(std::cos(phase2)*r1 + std::sin(phase2)*s1);
        mode.b[2]=amp*(std::cos(phase2)*r2 + std::sin(phase2)*s2);
        mode.norm = 0.5*((mode.a[0]*mode.a[0]+mode.a[1]*mode.a[1]+mode.a[2]*mode.a[2])+
                         (mode.b[0]*mode.b[0]+mode.b[1]*mode.b[1]+mode.b[2]*mode.b[2]));
        turb_modes.push_back(mode);
    }
    turb_coeff_init.assign(turb_modes.size(), 1.0);
    turb_coeff_curr.assign(turb_modes.size(), 0.0);
    turb_force_alpha = 1.0;
}

void sim_manager::isotropic_turbulence_velocity(double x, double y, double z,
        double &u, double &v, double &w) const {
    u = v = w = 0;
    if(!turb_iso_enable || turb_modes.empty()) return;
    for(size_t i=0;i<turb_modes.size();i++) {
        double bu, bv, bw;
        isotropic_turbulence_basis(i, x, y, z, bu, bv, bw);
        u += turb_coeff_init[i]*bu;
        v += turb_coeff_init[i]*bv;
        w += turb_coeff_init[i]*bw;
    }
}

void sim_manager::isotropic_turbulence_basis(size_t mode_id, double x, double y, double z,
        double &u, double &v, double &w) const {
    u = v = w = 0;
    if(mode_id >= turb_modes.size()) return;
    const turb_mode &m = turb_modes[mode_id];
    double phase = 2.0*M_PI*(m.k[0]*(x-ax)/lx + m.k[1]*(y-ay)/ly + m.k[2]*(z-az)/lz);
    double c = std::cos(phase);
    double s = std::sin(phase);
    u = m.a[0]*c + m.b[0]*s;
    v = m.a[1]*c + m.b[1]*s;
    w = m.a[2]*c + m.b[2]*s;
}

void sim_manager::isotropic_turbulence_force(double x, double y, double z,
        double &fx, double &fy, double &fz) const {
    if(!turb_iso_enable || turb_modes.empty()) return;
    for(size_t i=0;i<turb_modes.size();i++) {
        double bu, bv, bw;
        isotropic_turbulence_basis(i, x, y, z, bu, bv, bw);
        // force produces u_new = alpha*u_low + u_high over one step
        fx += turb_force_alpha*turb_coeff_curr[i]*bu;
        fy += turb_force_alpha*turb_coeff_curr[i]*bv;
        fz += turb_force_alpha*turb_coeff_curr[i]*bw;
    }
}

void sim_manager::update_isotropic_turbulence_state(double kinetic_energy, const std::vector<double> &proj_dot, double dt, double t){
    if(!turb_iso_enable) return;
    if(turb_modes.empty()) return;
    if(dt <= 0) return;
    if(proj_dot.size() != turb_modes.size()) return;

    double energy_lower = 0.0;
    for(size_t i=0;i<turb_modes.size();i++) {
        double norm = std::max(turb_modes[i].norm, 1e-14);
        turb_coeff_curr[i] = proj_dot[i] / norm;
        energy_lower += 0.5*turb_coeff_curr[i]*turb_coeff_curr[i]*norm;
    }
    double energy_upper = kinetic_energy - energy_lower;
    double alpha2 = (turb_target_energy - energy_upper) / std::max(energy_lower, 1e-12);
    if(alpha2 < 0) alpha2 = 0;
    double alpha = std::sqrt(alpha2);
    turb_force_alpha = (alpha-1.0)/dt;
    turb_last_update_t = t;
}

void sim_manager::add_obj(object *o) {
	// new array, copy old stuff over
	object **oarr = new object*[n_obj+1];
	for (int i = 0; i < n_obj; i++) oarr[i] = objs[i];

	// delete old array (if it exists)
	if (objs != NULL) delete[] objs;

	// assign new array to our pointer, set new val
	// and increment number of objects
	objs = oarr;
	objs[n_obj++] = o;
}

double sim_manager::phi(int obj_index, const double (&xi)[3]) const{
    if(obj_index<0 || obj_index >= n_obj) return nan("");
	return objs[obj_index]->phi(xi);
}

void sim_manager::velocity(double x,double y,double z,
	double &u,double &v,double &w) {
    u=0; v=0; w=0;
    const double coords[3] = {x,y,z};
    double rx[3] = {0,0,0};

    double * vpp = spars->vel_profile;
    for(int i=0;i<spars->vel_prof_num;i++){
        int vtype = static_cast<int> (vpp[2*i]);

        double amp = vpp[2*i+1];
        if(vtype == 0) continue;
        else if (vtype == 1) u+= amp;
        else if (vtype == 2) v+= amp;
        else if (vtype == 3) w+= amp;
        else {

        }
    }

    for(int i=0;i<n_obj;i++) {
        objs[i]->rm0(coords, rx);
        if(objs[i]->phi(rx) < 0.) {
            objs[i]->velocity(rx, u, v, w);
        }
    }

   if(turb_iso_enable){
        double tu, tv, tw;
        isotropic_turbulence_velocity(x, y, z, tu, tv, tw);
        u += tu;
        v += tv;
        w += tw;
    }
    
    double tmp = sqrt(u*u+v*v+w*w);
    if(tmp>vmax) vmax = tmp;
}

void sim_manager::solid_acceleration(const int obj_id, const double (&rval)[3], double x, double y, double z, double t,	double &fx,double &fy,double &fz){
//void sim_manager::solid_acceleration(const int obj_id, const double (&rval)[3], double x, double y, double z, double t,	double &fx,double &fy,double &fz, const double vx, const double vy, const double vz){
    fx=fy=fz=0;
}

// FLUID CONVERGENCE
void sim_ftest::velocity(double x,double y,double z,
	double &u,double &v,double &w) {
	u =     - sin(2*M_PI*x) * cos(2*M_PI*y) * cos(2*M_PI*z);
	v = 0.5 * cos(2*M_PI*x) * sin(2*M_PI*y) * cos(2*M_PI*z);
	w = 0.5 * cos(2*M_PI*x) * cos(2*M_PI*y) * sin(2*M_PI*z);
    double tmp = sqrt(u*u+v*v+w*w);
    if(tmp>vmax) vmax = tmp;
}

void sim_ftest::pressure(double x,double y,double z,double &p) {
	p = 3*M_PI*fm.mu*cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
}

void sim_ftest::exact(double x,double y,double z,double t,
	double h,double &u,double &v,double &w,double &p) {

	// get exact vels at cell center
	velocity(x,y,z,u,v,w);
	u *= cos(2*M_PI*t);
	v *= cos(2*M_PI*t);
	w *= cos(2*M_PI*t);

	// shift to corner for pressure
	x -= 0.5*h;
	y -= 0.5*h;
	z -= 0.5*h;
	pressure(x,y,z,p);
	p *= (cos(2*M_PI*t) - sin(2*M_PI*t)/(6*M_PI*fm.mu));
}

void sim_ftest::fluid_acceleration(double x,double y,double z,double t,
	double &fx,double &fy,double &fz) {
	fx = (12*M_PI*cos(2*M_PI*y)*cos(2*M_PI*z)*(-6*M_PI*fm.mu*cos(2*M_PI*t) +
            sin(2*M_PI*t))*sin(2*M_PI*x) + M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (2+cos(4*M_PI*y) + cos(4*M_PI*z))*sin(4*M_PI*x))/4;
    fy = -(M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (-1+cos(4*M_PI*x) -2*cos(4*M_PI*z))*sin(4*M_PI*y))/8;
    fz = -(M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (-1+cos(4*M_PI*x) -2*cos(4*M_PI*y))*sin(4*M_PI*z))/8;
}

// SOLID CONVERGENCE / SHEAR WAVE
sim_stest::sim_stest(const sim_params *spars) : sim_manager(spars),
	a(0.01),k(2*M_PI),w(2*M_PI) {
	shear_cube *cube = new shear_cube(a,k);
	add_obj(cube);
}

void sim_stest::exact(double x,double y,double z,double t,
	double h,double &u,double &v,double &w,double &p) {
	u = this->w*a*sin(k*z - w*t);
	v = this->w*a*sin(k*z - w*t);
	w = 0;

	x -= 0.5*h;
	y -= 0.5*h;
	z -= 0.5*h;
	p = -sm_array[0].G * a*a*k * (-k*cos(2*(k*z-this->w*t)) + cos(2*this->w*t - k)*sin(k))/3;
}

void sim_stest::pressure(double x,double y,double z,double &p) {
	p = -sm_array[0].G * a*a*k * (-k*cos(2*k*z) + cos(-k)*sin(k))/3;
}

void sim_stest::velocity(double x,double y,double z,
	double &u,double &v,double &w) {
	u = this->w*a*sin(k*z);
	v = this->w*a*sin(k*z);
	w = 0;
    double tmp = sqrt(u*u+v*v+w*w);
    if(tmp>vmax) vmax = tmp;
}

void sim_objects::solid_acceleration(const int obj_id, const double (&rval)[3], double x, double y, double z, double t,	double &fx,double &fy,double &fz){
        double krepx, krepy, krepz;
        krepx = wall_acc_mult*avg_velx[obj_id]*avg_velx[obj_id];
        krepy = wall_acc_mult*avg_vely[obj_id]*avg_vely[obj_id];
        krepz = wall_acc_mult*avg_velz[obj_id]*avg_velz[obj_id];

        if(krepx<min_wall_acc) krepx = min_wall_acc;
        if(krepy<min_wall_acc) krepy = min_wall_acc;
        if(krepz<min_wall_acc) krepz = min_wall_acc;

        // We use fx,fy,fz, but they are really accelerations
        fx=fy=fz=0;
        double phiv = phi(obj_id, rval);
        // wall forces in x direction
        if(walls[0]){
            double dist = x-wall_pos[0];
            if(dist<=wall_dist){
                fx += delta_func_wall(dist*winv) * krepx *wall_trans_func(phiv,dx,dxsp);
            }
        }
        if(walls[1]){
            double dist = wall_pos[1]-x;
            if(dist<=wall_dist){
                fx -= delta_func_wall(dist*winv) * krepx *wall_trans_func(phiv,dx,dxsp);
            }
        }

        // wall forces in y direction
        if(walls[2]){
            double dist = y-wall_pos[2];
            if(dist<=wall_dist){
                fy += delta_func_wall(dist*winv) * krepy *wall_trans_func(phiv,dy,dysp);
            }
        }
        if(walls[3]){
            double dist = wall_pos[3]-y;
            if(dist<=wall_dist){
                fy -= delta_func_wall(dist*winv) * krepy *wall_trans_func(phiv,dy,dysp);
            }
        }

        // wall forces in z direction
        if(walls[4]){
            double dist = z-wall_pos[4];
            if(dist<=wall_dist) {
                fz += delta_func_wall(dist*winv) * krepz *wall_trans_func(phiv,dz,dzsp);
            }
#if 0
            // HERE TO HARD CODE A WALL AT 1/2 slope, i.e. z = 0.5*x
            double dist = z-0.5*x;
            // normal distance from a point to the plane that has a 1/2 slope
            double norm_dist = dist * 0.894;
            if(norm_dist<=wall_dist){
                double tmp_force = delta_func_wall(dist*winv) * krepz *wall_trans_func(phiv,dz,dzsp);
                fz += tmp_force * 0.8;
                fx -= tmp_force * 0.4;
            }
#endif
        }
        if(walls[5]){
            double dist = wall_pos[5]-z;
            if(dist<=wall_dist){
                fz -= delta_func_wall(dist*winv) * krepz *wall_trans_func(phiv,dz,dzsp);
            }
        }
        double f_norm = sqrt(fx*fx + fy*fy + fz*fz);
        if(f_norm >= K_stiff ){
            double normfac = K_stiff/f_norm;
            fx *= normfac;
            fy *= normfac;
            fz *= normfac;
        }
        // Sometimes the anchors edge need to be specified by 2 distances
        // e.g. a cylindrical anchor region, R, L
        double ldx=0, ldy=0, ldz=0, dtoe1=2*eps, dtoe2=2*eps;
        objs[obj_id]->passive_acc_deformation(rval, eps, x, y, z, t, ldx, ldy, ldz, dtoe1, dtoe2);
        double K = K_stiff*heaviside(dtoe1)*heaviside(dtoe2);
        fx-=K*ldx;
        fy-=K*ldy;
        fz-=K*ldz;

        // Add gravitational acceleration
        fz += gravity*heaviside(phiv);
}

// FLUID CONVERGENCE
void object_mills::velocity(double x,double y,double z,
	double &u,double &v,double &w) {
	u =     - sin(2*M_PI*x) * cos(2*M_PI*y) * cos(2*M_PI*z);
	v = 0.5 * cos(2*M_PI*x) * sin(2*M_PI*y) * cos(2*M_PI*z);
	w = 0.5 * cos(2*M_PI*x) * cos(2*M_PI*y) * sin(2*M_PI*z);
    double tmp = sqrt(u*u+v*v+w*w);
    if(tmp>vmax) vmax = tmp;
}

void object_mills::pressure(double x,double y,double z,double &p) {
	p = 3*M_PI*fm.mu*cos(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
}

void object_mills::exact(double x,double y,double z,double t,
	double h,double &u,double &v,double &w,double &p) {

	// get exact vels at cell center
	velocity(x,y,z,u,v,w);
	u *= cos(2*M_PI*t);
	v *= cos(2*M_PI*t);
	w *= cos(2*M_PI*t);

	// shift to corner for pressure
	x -= 0.5*h;
	y -= 0.5*h;
	z -= 0.5*h;
	pressure(x,y,z,p);
	p *= (cos(2*M_PI*t) - sin(2*M_PI*t)/(6*M_PI*fm.mu));
}

void object_mills::fluid_acceleration(double x,double y,double z,double t,
	double &fx,double &fy,double &fz) {
	fx = (12*M_PI*cos(2*M_PI*y)*cos(2*M_PI*z)*(-6*M_PI*fm.mu*cos(2*M_PI*t) +
            sin(2*M_PI*t))*sin(2*M_PI*x) + M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (2+cos(4*M_PI*y) + cos(4*M_PI*z))*sin(4*M_PI*x))/4;
    fy = -(M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (-1+cos(4*M_PI*x) -2*cos(4*M_PI*z))*sin(4*M_PI*y))/8;
    fz = -(M_PI*cos(2*M_PI*t)*cos(2*M_PI*t) *
            (-1+cos(4*M_PI*x) -2*cos(4*M_PI*y))*sin(4*M_PI*z))/8;
}

