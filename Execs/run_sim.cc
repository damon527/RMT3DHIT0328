#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

#include "sim_manager.hh"
#include "geometry.hh"
#include "stepper.hh"
#include "fluid_3d.hh"
#include "sim_params.hh"

sim_manager * sim_manager_alloc(sim_params *spars){
    const char * type = spars->sim_type;
    if (strcmp(type,"ftest") == 0) {
        return new sim_ftest(spars);
    } else if (strcmp(type,"stest") == 0) {
        return new sim_stest(spars);
    } else if (strcmp(type,"objects") == 0) {
        return new sim_objects(spars);
    } else if (strcmp(type,"object_mills") == 0) {
        return new object_mills(spars);
    }
    else{
        printf("Simulation types can be:\n(1) ftest\n(2) stest\n(3) objects\n");
		printf("Provided: |%s|\n",type);
        printf("Modify your .cfg file and try again. A complete set of parameters are printed for you in %s/recover.cfg\n", spars->dirname);
        spars->write_params(spars->dirname);
        MPI_Abort(world, 1);
        return NULL;
    }
}

/**
 * Driver function for the fluid_3d class, solving incompressible
 * Navier-Stokes with Chorin projection method as implemented in
 * Yu inkjet paper
 */
int main(int argc, char* argv[]) {

	// Initiate MPI, get rank and the number of processes
	int rank,procs;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(world,&rank);
	MPI_Comm_size(world,&procs);

	// ########### input validation ############
	if(argc<2) {
		printf("Usage: mpirun -np N [-mpiargs...] run_sim sample_sim.cfg\n");
		MPI_Finalize();
		exit(0);
	}

    // Print a time and date stamp
    char * system_cmd = NULL;
    if(rank==0) {
        system_cmd = new char[256];
        sprintf(system_cmd, "echo \"Simulation started at `date`\"");
        int err = system(system_cmd);
        if(err) printf("Oops, system call data didn't work\n");
    }


	sim_params spars(argv[1]);
	geometry gm(spars.sys_size[0],spars.sys_size[1],spars.sys_size[2],spars.x_prd,spars.y_prd,spars.z_prd);
	// initialize sim_type and geometry, get CFL max dt
	sim_manager *mgmt = sim_manager_alloc(&spars);
	  const int configured_n_obj = spars.n_obj;
    const bool do_hit_cold_start = spars.hit_cold_start_enable && spars.turb_iso_enable && (configured_n_obj>0);
    const bool use_recovered_background = do_hit_cold_start && (spars.chk_num>0);
    if(do_hit_cold_start && !use_recovered_background) {
        // Single-phase warmup stage.
        mgmt->n_obj = 0;
    }
    if(rank==0) mkdir(spars.dirname,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

	if (rank == 0) {
		printf("\n#**************** %2d procs *****************\n"
			"# Simulation type: %s\n"
			"# initialized using %d x %d x %d decomposition\n"
			"# global grid size %d x %d x %d\n"
			"# global domain size (%3.2f, %3.2f)x(%3.2f, %3.2f)x(%3.2f, %3.2f)\n"
			"# periodicity %d x %d x %d\n"
			"# output directory: %s\n#\n"
			,procs,
			spars.sim_type,
			gm.mp,gm.np,gm.op,
			spars.sys_size[0],
			spars.sys_size[1],
			spars.sys_size[2],
            spars.ax, spars.bx,
            spars.ay, spars.by,
            spars.az, spars.bz,
            spars.x_prd,spars.y_prd,spars.z_prd,
            spars.dirname);
        printf("# Simulating %d objects\n", spars.n_obj);
        for(int i=0;i<spars.n_obj;i++) mgmt->objs[i]->print_self();
    }

	// now create a fluid_3d simulator
	fluid_3d f3d(&gm, mgmt, &spars);
    // Part of initialization also sets time scale, espeically when there's non-trivial initial velocities
    int init_err = f3d.initialize();
    f3d.setup_hit_background();

	// check if user specified time spacing, combine with
	// frames for stepper object
	if(mgmt->dt_reg<spars.dt){
        spars.dt = mgmt->dt_reg;
    }

f3d.init_iter(init_err);

    if(do_hit_cold_start) {
        if(rank==0) {
            if(use_recovered_background) {
                printf("# HIT cold-start: recovered background detected, skip single-phase warmup.\n");
            } else {
                printf("# HIT cold-start: warm up single-phase HIT for %.6g seconds.\n", spars.hit_cold_start_time);
            }
        }

        if(!use_recovered_background) {
            const int warm_steps = (spars.hit_cold_start_time>0. && spars.dt>0.) ? static_cast<int>(ceil(spars.hit_cold_start_time/spars.dt)) : 0;
            if(rank==0) printf("# HIT cold-start: warmup steps = %d\n", warm_steps);
            for(int i=0;i<warm_steps;i++) f3d.step_forward(spars.debug_flag);

            if(spars.hit_write_background_chk) {
                int bg_chk = spars.hit_background_chk_num;
                if(bg_chk<0) bg_chk = 0;
                char bg_chk_dir[256];
                sprintf(bg_chk_dir, "%s/chk.%05d", spars.dirname, bg_chk);
                if(rank==0) mkdir(bg_chk_dir,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
                f3d.write_chk_pt(warm_steps, bg_chk, bg_chk_dir);
                if(rank==0) printf("# HIT cold-start: wrote background checkpoint %s\n", bg_chk_dir);
            }
        }

        // Insert particles and optionally keep/stop HIT forcing.
        f3d.set_active_object_count(configured_n_obj);
        f3d.insert_particles_from_current_objects();
        f3d.set_hit_forcing_enabled(spars.hit_keep_forcing_after_insert);
        if(rank==0) {
            printf("# HIT cold-start: particle insertion done, forcing after insert = %s\n",
                   spars.hit_keep_forcing_after_insert?"ON":"OFF");
        }

        const int pre_steps = (spars.hit_pre_iter_steps>0)?spars.hit_pre_iter_steps:0;
        if(rank==0) printf("# HIT cold-start: pre-iterations after insertion = %d\n", pre_steps);
        for(int i=0;i<pre_steps;i++) f3d.step_forward(spars.debug_flag);

        // Preprocessing complete: reset counters for production run.
        f3d.nt = 0;
        f3d.time = 0;
        spars.chk_step = 0;
        spars.chk_num = 0;
        spars.set_current_time(0);
    }

	stepper sim(spars,f3d);
	f3d.dt = sim.dt;

    if(rank == 0){

        mgmt->fm.print_self();
        for(int i=0; i<mgmt->n_obj; i++){
            mgmt->sm_array[i].print_self();
        }
        mgmt->print_self();
        printf("# Simulating to T = %0.3g with dt = %.10e (double check %.10e)\n", spars.T, spars.dt, sim.dt);

        if (spars.dump_code){
            printf("# Start at %d frames, output %d frames of %d steps.\n"
                   "# Last frame has additional %d steps.\n",
                sim.init_frame, sim.frames, sim.steps, sim.last_steps);
                if(spars.dump_code & 1) { puts("# Output tracers."); }
                if(spars.dump_code & 2) { puts("# Output contours."); }
                if(spars.dump_code & 4) {
                    puts("# Output slices.");
                    printf("# Slice flag %d, format flag %d, slice in %d dir at %d\n"
                        , spars.out_flag, spars.data_fmt, spars.output_dim, spars.output_ind);
                    for(int i=0;i<write_params::numf;i++) {
                        if((1<<i) & spars.out_flag) {
                            printf("#   Output field %s\n", write_params::default_names[i]);
                        }
                    }
                }
                if(spars.dump_code & 8) { puts("# Output checkpoint files."); }
                if(spars.dump_code & 16) { puts("# Output diagonastics file."); }
                if(spars.dump_code & 32) { puts("# Output macro stats file."); }

        }
	}

	// run simulation
	sim.run(f3d);
	//f3d.compute_stress();
	// print error if want it
	if (mgmt->has_exact) {
		double err[4];
		f3d.error(err);
		if (rank == 0) {
			//double t[2];
			//int u[2];
			//for(int i=0;i<2;i++) u[i] = f3d.watch.convert(t[i], i);
			printf("\n\n#T dh dt erru_sq errv_sq errw_sq errp_sq\n");
			printf("%16.14g\t%13.10f\t%13.10f\t%13.10f\t%13.10f\t%13.10f\t%13.10f\n\n\n",
				f3d.time, 1.0/spars.sys_size[0], spars.dt, err[0],err[1],err[2],err[3]);
			//printf("\nOld extrapolation took %11.8f %c\n"
			//"New extrapolation took %11.8f %c\n",
			//	t[0], u[0], t[1], u[1]);
		}
	}
    if(rank==0){
        f3d.watch.report();
        f3d.expper.watch.report();
    }

	gm.free();
	f3d.free();

    if(rank==0) {
        sprintf(system_cmd, "echo \"Simulation ended at `date`\"");
        int err = system(system_cmd);
        if(err) printf("Oops, system call data didn't work\n");
        delete [] system_cmd;
    }
    delete mgmt;
	MPI_Finalize();
}
