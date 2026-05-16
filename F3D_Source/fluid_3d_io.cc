#include "fluid_3d.hh"

void fluid_3d::write_slice(write_params wp,const char* filename) {
	const int dinfo[5]={2,2,1,0,0};
	int &dim=wp.dim;

	// get global dimensions, setting slice coord to 1
	int global_dim[3] = {m, n, o};
	// if there are four items in the field, it's fluid field, hack for now
	if(wp.corner_field()) {
		global_dim[0] = fem_grid.m;
		global_dim[1] = fem_grid.n;
		global_dim[2] = fem_grid.o;
	}
	global_dim[dim] = 1;

	// find out the other two dimension to iterate through
	int iter_dim1 = dinfo[dim+2], iter_dim2=dinfo[dim];

	int gslice_len = global_dim[0] * global_dim[1] * global_dim[2];
	// declare array, instantiate if master, and pull out slice

	double *u_global = NULL;
	if (rank == 0) {
		u_global = new double [gslice_len];
		for (int i=0; i<gslice_len; i++) *(u_global+i) = 0;
	}
	slice(wp, u_global);

	// if we're not the master, we're done. wait up and then quit.
	if (rank != 0) {
		MPI_Barrier(grid->cart);
		return;
	}

	// if we're the master, open file...
	FILE *fh = p_safe_fopen(filename, "w");

	// march through writing file (starting with largest y)
	if(wp.format==0) save_matrix(fh,u_global,global_dim[iter_dim1],global_dim[iter_dim2]);
	else {
		double *coord1 = gx;
		double *coord2 = gy;
		if(dim==0) { coord1 = gy; coord2 = gz; }
		else if(dim==1) { coord2 = gz; }
		if(wp.format==1) save_text(fh,u_global,global_dim[iter_dim1],global_dim[iter_dim2],coord1,coord2);
		else save_gnuplot(fh,u_global,global_dim[iter_dim1],global_dim[iter_dim2],coord1,coord2);
	}

	// close file and deallocate array
	fclose(fh);
	delete[] u_global;

	// join barrier, then we're done
	MPI_Barrier(grid->cart);
}

void fluid_3d::save_matrix(FILE *fh,double *u_global,int d1,int d2) {
	for (int j=0;j<d2;j++) {
		for (int i = 0; i < d1; i++) fprintf(fh,"%17.14g ", u_global[i + j*d1]);
		fputc('\n', fh);
	}
}

void fluid_3d::save_text(FILE *fh,double *u_global,int d1,int d2,double *a1,double *a2) {
	for (int j=0;j<d2;j++) {
		for (int i = 0; i < d1; i++)
			fprintf(fh,"%17.14f %17.14f %17.14f\n", a1[i], a2[j], u_global[i + j*d1]);
		fputc('\n', fh);
	}
}

void fluid_3d::save_gnuplot(FILE *fh,double *u_global,int d1,int d2,double *a1,double *a2) {
	float *buf=new float[d1+1],*bp=buf+1;*buf=static_cast<float>(m);
	int i,j;

	// Output the size and x coordinates
	for(i=0;i<d1;i++) bp[i]=static_cast<float>(a1[i]);
	fwrite(buf,sizeof(float),d1+1,fh);

	// Output the y coordinate and field values
	double *up=u_global;
	for(j=0;j<d2;j++) {
		*buf=static_cast<float>(a2[j]);
		for(i=0;i<d1;i++) bp[i]=static_cast<float>(*(up++));
		fwrite(buf,sizeof(float),d1+1,fh);
	}
	delete [] buf;
}

/**
* Grab a slice of nodes at constant (global) x-, y-, or z-index
* position. Only on master node, will write into array at address
* g_val in the form u(xi,yj) -> u_global[xi + m*yj] if slicing z.
* \param[in] wp the parameters for this output.
* \param[in] g_val the pointer to write out to.
*/
void fluid_3d::slice(write_params wp, double *g_val) {
	int &dim = wp.dim, &point = wp.point;

	bool self_involved=false;
	double * send_recv_buf = cbuf.buf;
	// Put needed local and global into array of 3 elements
	// in the order of x, y, z
	const int global_dim [3] = {m,n,o};
	const int p_global_dim [3] = {fem_grid.m,fem_grid.n,fem_grid.o};
	const int local_dim [3] = {sm,sn,so};
	const int p_local_dim [3] = {fem_grid.sm,fem_grid.sn,fem_grid.so};
	const int num_proc [3] ={grid->mp,grid->np,grid->op};
	const int proc_id[3] = {grid->ip,grid->jp,grid->kp};
	const int lowest [3] = {ai,aj,ak};
	const int highest[3] = {bi,bj,bk};

	// find out the other two dimension to iterate through
	int tmp1 = (dim + 1) % 3, tmp2 = (dim + 2) % 3;

	// iter_dim1 is the smaller of the two dimensions
	// suppose we are doing z slice, iter_dim1=x, iter_dim2 = y
	// supoose we are doing y slice, iter_dim1=x, iter_dim2 = z
	const int iter_dim1 = (tmp1<tmp2)?tmp1:tmp2;
	const int iter_dim2 = (tmp1<tmp2)?tmp2:tmp1;
	// make array of sender, length, dimension in global array
	int nprocs1 = num_proc[iter_dim1];
	int nprocs2 = num_proc[iter_dim2];
	// needed for data assembly
	int global_len = global_dim[iter_dim1];
	int local_dim1 = local_dim[iter_dim1];
	int local_dim2 = local_dim[iter_dim2];

	if(wp.corner_field()) {
		global_len = p_global_dim[iter_dim1];
		local_dim1 = p_local_dim[iter_dim1];
		local_dim2 = p_local_dim[iter_dim2];
	}

	// array and field labels
	int tag_n_len [5];
	const int STAG = 0;
	const int SLEN = 1;
	const int LDIM1 = 2;
	const int LDIM2 = 3;
	const int GDIM1 = 4;
	const int SENDR = 5;

	// find the processor index containing slice along perpendicular dimension
	int slice_ind = -1;
	if(point>=lowest[dim] && point<highest[dim]) slice_ind=proc_id[dim];

    // The following 4 lines is so that if the master node doesn't belong to the slice
    // it will still know which index the slice is at
	double tmp_ind_local = (slice_ind>=0)?((double) slice_ind / double (nprocs1*nprocs2)):0;
	double tmp_ind_master;
	MPI_Reduce(&tmp_ind_local, &tmp_ind_master, 1, MPI_DOUBLE, MPI_SUM, 0, grid->cart);
	if(rank==0 && slice_ind==-1) slice_ind=(int) tmp_ind_master;

	// just random tag to send tag
	const int TAGTAG = 2;
	int slice_len, slice_tag;

	// if this is a processor that contains the slice in question...
    if(wp.o_type == 18) {
        fill_pp_rhs();
    }
	if (proc_id[dim] == slice_ind) {
		// load the info on to the writing array
		slice_len = copy_slice_to_buf(wp,point-lowest[dim]);
		// the follow algebraic gymnastic to figure out the tags
		// tag is the starting position each local slice
		// but the position is computed with respect to the global slice
		// for a z slice, my tag is ai+aj*m
		// for a y slice, my tag is ai+ak*m
		// for a x slice, my tag is aj+ak*n
		slice_tag = lowest[iter_dim1] + lowest[iter_dim2] * global_len;
		// First send the slice tags to master node, with the tag -1
		tag_n_len[STAG] = slice_tag;
		tag_n_len[SLEN] = slice_len;
		tag_n_len[LDIM1] = local_dim1;
		tag_n_len[LDIM2] = local_dim2;
		tag_n_len[GDIM1] = global_len;

		if (rank==0) {
			// since master node doesn't send stuff to itself
			// we manually increment the buffer pointer, ready to receive
			send_recv_buf += slice_len;
			self_involved=true;
		}
		else {
			MPI_Send(tag_n_len,5,MPI_INT,0,TAGTAG,grid->cart);
			// then send the slice data to master node with slice tags
			MPI_Isend(send_recv_buf, slice_len, MPI_DOUBLE, 0, slice_tag, grid->cart, reqs);
			MPI_Waitall(1,reqs,stats);
		}
	}
	// now that all's sent, we collect on the master
	if (rank == 0) {

		// define array for cartesian processor coords
		int coords[3] = {0,0,0};
		coords[dim] = slice_ind;

		int **send_list = new int* [nprocs1*nprocs2];
		// now for each processor
		for (int i = 0; i < nprocs1*nprocs2; i++) {
			// make array of size six because we also store the order
			// in which these slice messages were received, i.e. SENDR
			send_list[i] = new int[6];
			for (int j = 0; j<6; j++) send_list[i][j] = 0;
		}
		int ** sublist = send_list;
		if(self_involved){
			// we in fact did not send from master to master
			// manually load data
			for (int i=0;i<5;i++) send_list[0][i] = tag_n_len[i];
			send_list[0][5] = 0;
			sublist++;
		}

		// iterate through all the processors that intersect
		// the slice, and recv tag and length info
		int sender, comm_co=0;
		for(int i = 0; i < nprocs1; i++) {
			for(int j = 0; j < nprocs2; j++) {

				// record the processor coords and get rank
				coords[iter_dim1] = i;
				coords[iter_dim2] = j;
				MPI_Cart_rank(grid->cart,coords,&sender);
				if (sender==0) continue;

				// get pointer to list for storing info
				*(*sublist+SENDR) = sender;
				// fill in rest of list with info from sender
				MPI_Recv(*sublist,5,MPI_INT,sender,TAGTAG,grid->cart,MPI_STATUS_IGNORE);
				// printf("from sender %d, sublist : { ", sender); for (int k=0;k<5;k++) printf("%d ", *(*sublist+k)); printf(" }\n");
				// includes slice data tag and length of slice data
				int pos = *(*sublist+STAG);
				int len = *(*sublist+SLEN);
				// now pull slice data into buffer
				MPI_Irecv(send_recv_buf,len,MPI_DOUBLE,sender,pos,grid->cart,reqs+comm_co);
				send_recv_buf += len;
				sublist++;comm_co++;
			}
		}
		// wait to finish everything, then copy out of buffer
		MPI_Waitall(comm_co,reqs,stats);
		copy_slice_from_buf(send_list,nprocs1*nprocs2,g_val);

		// delete the new array we made
		for (int i = 0; i < nprocs1*nprocs2; i++) if (send_list[i]!=NULL) delete [] send_list[i];
		if(send_list!=NULL) delete[] send_list;
	}
}

/** Copy slice at const i at local index ind into buffer to send.
* \param[in] wp the parameters for this output.
* \param[in] point the local coordinate of the slice to print. */
int fluid_3d::copy_slice_to_buf(write_params wp,int point) {
	int &dim = wp.dim;

	//printf("I'm processor %d, (%d, %d, %d), slicing at %d, point %d\n", rank, ai, aj, ak, dim, point);
	if(point<0 || (dim==0 && point>=sm) || (dim==1&&point>=sn) || (dim==2 && point>=so)){
		printf("copy_slice_to_buf(): Error! Index out of bound. Nothing has been done.\n");
		return 0;
	}

	// set limits for for loop -
	// start (change dir. in question to fixed point)
	int sp[3] = {0,0,0};

	// end position (dir. in question to fix_point + 1)
	int ep[3] = {sm, sn, so};
	// if writing out pressure we need to use fem_grid dimensions
	if(wp.corner_field()) {
		ep[0] = fem_grid.sm;
		ep[1] = fem_grid.sn;
		ep[2] = fem_grid.so;
	}

    sp[dim] = point;
	ep[dim] = point+1;
	double *b = cbuf.buf;
	for (int k = sp[2]; k<ep[2]; k++) for (int j = sp[1]; j < ep[1]; j++)
		for (int i = sp[0]; i < ep[0]; i++, b++) {
		// not outputing levelset
		int eid = i + sm4*j + smn4*k;
		int eidlow = eid-smn4;

        int obj_id=-1;
		if(wp.o_type<3) {
			switch(wp.o_type) {
				case 0: *b=u0[eid].vel[0];break;
				case 1: *b=u0[eid].vel[1];break;
				case 2: *b=u0[eid].vel[2];break;
            }
        } else if(wp.o_type==3) {
                // wx d(vz)/dy - d(vy)/dz
                double tmp_divu = 0.25*dysp*(u0[eid].vel[2] - u0[eid-sm4].vel[2] + u0[eid-1].vel[2] - u0[eid-1-sm4].vel[2]
				           + u0[eidlow].vel[2] - u0[eidlow-sm4].vel[2] + u0[eidlow-1].vel[2] - u0[eidlow-1-sm4].vel[2])
				- 0.25*dzsp*(u0[eid].vel[1] - u0[eidlow].vel[1] + u0[eid-1].vel[1] - u0[eidlow-1].vel[1]
					   + u0[eid-sm4].vel[1] - u0[eidlow-sm4].vel[1] + u0[eid-1-sm4].vel[1] - u0[eidlow-1-sm4].vel[1]);
                *b = tmp_divu;
        } else if(wp.o_type==4) {
                // wy d(vx)/dz - d(vz)/dx
                double tmp_divu = 0.25*dzsp*(u0[eid].vel[0] - u0[eidlow].vel[0] + u0[eid-1].vel[0] - u0[eidlow-1].vel[0]
					   + u0[eid-sm4].vel[0] - u0[eidlow-sm4].vel[0] + u0[eid-1-sm4].vel[0] - u0[eidlow-1-sm4].vel[0])
				- 0.25*dxsp*(u0[eid].vel[2] - u0[eid-1].vel[2] + u0[eid-sm4].vel[2] - u0[eid-sm4-1].vel[2]
					   + u0[eidlow].vel[2] - u0[eidlow-1].vel[2] + u0[eidlow-sm4].vel[2] - u0[eidlow-sm4-1].vel[2]);
                *b = tmp_divu;
        } else if(wp.o_type==5) {
                // wz d(vy)/dx - d(vx)/dy
                double tmp_divu = 0.25*dxsp*(u0[eid].vel[1] - u0[eid-1].vel[1] + u0[eid-sm4].vel[1] - u0[eid-sm4-1].vel[1]
					   + u0[eidlow].vel[1] - u0[eidlow-1].vel[1] + u0[eidlow-sm4].vel[1] - u0[eidlow-sm4-1].vel[1])
				- 0.25*dysp*(u0[eid].vel[0] - u0[eid-sm4].vel[0] + u0[eid-1].vel[0] - u0[eid-1-sm4].vel[0]
				           + u0[eidlow].vel[0] - u0[eidlow-sm4].vel[0] + u0[eidlow-1].vel[0] - u0[eidlow-1-sm4].vel[0]);

                *b = tmp_divu;
        } else if(wp.o_type==6) {
                double tmp_magu = 0;
                for(int ii=0;ii<3;ii++) {
                        tmp_magu += u0[eid].vel[ii] * u0[eid].vel[ii];
                }
                *b = sqrt(tmp_magu);
        } else if(wp.o_type==7){
                *b = u0[eid].p;
		} else if (wp.o_type == 18){
            *b = mg_pres->reg[0]->r0[i + fem_grid.sm * j + fem_grid.sm * fem_grid.sn * k];
        } else if (wp.o_type == 19) {
            double tmp_divu = dxsp * (u0[eid].fvel[1][0] - u0[eid].fvel[0][0]) + dysp * (u0[eid].fvel[3][1] - u0[eid].fvel[2][1]) + dzsp * (u0[eid].fvel[5][2] - u0[eid].fvel[4][2]);
            *b = tmp_divu;
        } else if (wp.o_type == 20) {
#if defined(VAR_DEN)
            *b = lrho[eid+G0];
#else
            *b = 0.;
#endif
        } else if (wp.o_type == 21) {
#if defined (DEBUG)
            *b = lJ[eid+G0];
#else
            *b = 100.;
#endif
        } else {
			if(wp.obj_id==-1) {
				obj_id=min_phi_id(eid);
				if(obj_id==-1) {
					if(wp.o_type ==12){
						*b=0.;
					} else if(wp.o_type<12 || wp.o_type > 12){
						*b=100.;
					}
					continue;
				}
			} else obj_id=wp.obj_id;

			if(rm0[eid].oid() == obj_id) {
				switch(wp.o_type) {
					case 8: *b = rm0[eid].phi(mgmt);break;
					case 9: *b = rm0[eid].x[0];    break;
					case 10: *b = rm0[eid].x[1];   break;
					case 11: *b = rm0[eid].x[2];   break;
					case 12: *b = static_cast<double> (rm0[eid].oid()); break;
					case 13: *b = static_cast<double> (rm0[eid].lid()); break;
					case 14: *b = rm0[eid].phi(mgmt,1);      break;
					case 15: *b = rm0[eid].xpred[0];        break;
					case 16: *b = rm0[eid].xpred[1];        break;
					case 17: *b = rm0[eid].xpred[2];        break;
				}
			} else {
				bool found=false;
				/*
				if(eid==1017723 || eid==1017725 || eid==1017926){
					printf("IO: Extrapolate maps at %d : %d\n", eid, extraps.n0[eid]);
				}
				*/
				for(int n=0;n<extraps.n0[eid];n++){
					ref_map &p = extraps.f0[eid][n];
					/*
						if(eid==1017723 || eid==1017725 || eid==1017926){
								printf("IO: map %d obj id %d layer id %d\n", n, p.oid(), p.lid());
						}
					*/
					if(p.oid() == obj_id) {
						switch(wp.o_type) {
							case 8: *b = p.phi(mgmt);break;
							case 9: *b = p.x[0];break;
							case 10: *b = p.x[1];break;
							case 11: *b = p.x[2];break;
                            // we are not printout object id for extrapolated values
                            // since object id is only used for checkpoint purposes
                            // and only primary grid refmap is needed for that
                            case 12: *b = 0.; break;
                            case 13: *b = static_cast<double> (p.lid()); break;
                            case 14: *b = p.phi(mgmt,1);break;
                            case 15: *b = p.xpred[0];  break;
                            case 16: *b = p.xpred[1];  break;
                            case 17: *b = p.xpred[2];  break;
						}
						found=true;break;
					}
				}
				if(!found) {
                    if(wp.o_type==12) *b=0.;
                    else *b=100.;
                }
			}
		}

	}

	// calculate and return length of slice copied to buffer
	int slice_len = 1;
	for (int i = 0; i < 3; i++) slice_len *= ep[i] - sp[i];
	return slice_len;
}

/** Computes the object ID for the minimum phi value at a gridpoint. If no phi
 * value is available, it returns -1.
 * \param[in] eid the gridpoint to consider.
 * return The object ID corresponding to the minimum phi value. */
int fluid_3d::min_phi_id(int eid) {
	int lo=rm0[eid].oid(),obj_id=-1;
	double mphi=1e30,tphi;
	if(lo>0) {obj_id=lo;mphi=rm0[eid].phi(mgmt);}
	for(int n=0;n<extraps.n0[eid];n++) {
		ref_map &p=extraps.f0[eid][n];
		tphi=p.phi(mgmt);
		if(tphi<mphi) {
			mphi=tphi;
			obj_id=p.oid();
		}
	}
	return obj_id;
}

/** Computes the minimum phi value a gridpoint. If no phi
 * value is available, it returns 100.
 * \param[in] eid the gridpoint to consider.
 * return The minimum phi value. */
double fluid_3d::min_phi(int eid) {
	double mphi=100,tphi;
	//printf("min_phi: eid %d\n", eid);
	if(rm0[eid].oid()>0) {
        //printf("Output: %d, x here %g\n", eid, rm0[eid].x[0]);
        mphi=rm0[eid].phi(mgmt);
    }
	for(int n=0;n<extraps.n0[eid];n++) {
		ref_map &p=extraps.f0[eid][n];
		tphi=p.phi(mgmt);
		if(tphi<mphi) mphi=tphi;
	}
	return mphi;
}

/** Computes the phi value a gridpoint corresponding to a specific object. If no
 * phi value is available, it returns 100.
 * \param[in] eid the gridpoint to consider.
 * \param[in] cnum the object to consider.
 * return The minimum phi value. */
double fluid_3d::phi(int eid,int cnum) {
	if(rm0[eid].oid()==cnum) return rm0[eid].phi(mgmt);
	for(int n=0;n<extraps.n0[eid];n++) {
		ref_map &p=extraps.f0[eid][n];
		if(p.oid()==cnum) return p.phi(mgmt);
	}
	return 100;
}

/* ##################### DATA OUTPUT #####################*/
/**
 * copy slice from buffer to the global array.
 * \param[in] send_list is an double pointer, indexed by the order
 *     messaged arrived. For each message, send_list[n] includes
 *     starting position [0], ldim1[2], ldim2[3], and gdim1[4].
 *     These are the starting position of the message in global
 *     array, dimensions of the local slice, low and high
 *     (respectively), and global low dimension.
 * \param[in] g_val pointer to write out to.
 */
void fluid_3d::copy_slice_from_buf(int **send_list, int co, double *g_val) {

	// fields of parameter array
	const int STAG = 0;
	const int LDIM1 = 2;
	const int LDIM2 = 3;
	const int GDIM1 = 4;

	// info for for loop
	int slice_tag,dim1,dim2,gdim1;

	// go to buffer
	double *b = cbuf.buf;
	for (int n = 0; n < co; n++) {
		// read in fields from array
		slice_tag = send_list[n][STAG];
		dim1 = send_list[n][LDIM1];
		dim2 = send_list[n][LDIM2];
		gdim1 = send_list[n][GDIM1];

		// loop through, copying from buffer
		for (int j = 0; j < dim2; j++)
			for(int i = 0; i < dim1; i++, b++) {
				*(g_val+slice_tag+i+j*gdim1) = *b;
		}
	}
}

void fluid_3d::set_dump(int dim,int pnt) {
	dd = dim;
	ds = pnt;
}

void fluid_3d::display_stats() {
	out = new char[1024];
	if (rank == 0) sprintf(out,"multigrid statistics\n");
}

void fluid_3d::dump(int num,char *dir) {
/*		char fid[4] = {'u','v','w','p'};
		char *file = new char[256];
		for (int f = 0; f < 4; f++) {
			sprintf(file,"%s/%c%d",dir,fid[f],num);
			write_params<field> params(dd, ds, 4, f, u0, file);
			write_slice<field>(params);
		}
		if (trace) {
			sprintf(file,"%s/t%d",dir,num);
			tr->print(file);
		}
		delete[] file;*/
}

/** Sets up the array that gives the dimensions of the other processors, needed
 * for output and also for grid transfers. */
void fluid_3d::setup_output_dimensions() {

	if(rank==0) gather_sizes();
	else {

		// On other ranks, check that there will be enough space to
		// send a local slice of output. For processors that are
		// orthogonally aligned with the (0,0,0) processor, send
		// information about the grid dimension.
		if(grid->kp==0) {
			if(grid->jp==0) MPI_Send(&sm,1,MPI_INT,0,msg_trans_dims,grid->cart);
			if(grid->ip==0) MPI_Send(&sn,1,MPI_INT,0,msg_trans_dims,grid->cart);
		} else if(grid->ip==0&&grid->jp==0) MPI_Send(&so,1,MPI_INT,0,msg_trans_dims,grid->cart);
	}
}

/** A routine run by the master processor to gather information about the
 * dimensions of the other regions. The routine also calculates the maximum
 * dimensions in each direction, which is needed for allocating memory for the
 * contour plotter. */
void fluid_3d::gather_sizes() {
	int q[4],&i=*q,&j=q[1],&k=q[2],l;j=k=0;
	int &mp=grid->mp,&np=grid->np,&op=grid->op;
	osm=new int[mp+np+op+3];osn=osm+mp;oso=osn+np;max_sizes=oso+op;
	int &msm=*max_sizes,&msn=max_sizes[1],&mso=max_sizes[2];

	// Receive dimensions in the x direction
	msm=*osm=sm;
	for(i=1;i<mp;i++) {
		MPI_Cart_rank(grid->cart,q,q+3);
		MPI_Recv(osm+i,1,MPI_INT,q[3],msg_trans_dims,grid->cart,stats);
		if(osm[i]>msm) msm=osm[i];
	}

	// Receive dimensions in the y direction
	msn=*osn=sn;i=0;
	for(j=1;j<np;j++) {
		MPI_Cart_rank(grid->cart,q,q+3);
		MPI_Recv(osn+j,1,MPI_INT,q[3],msg_trans_dims,grid->cart,stats);
		if(osn[j]>msn) msn=osn[j];
	}

	// Receive dimensions in the z direction
	mso=*oso=so;j=0;
	for(k=1;k<op;k++) {
		MPI_Cart_rank(grid->cart,q,q+3);
		MPI_Recv(oso+k,1,MPI_INT,q[3],msg_trans_dims,grid->cart,stats);
		if(oso[k]>mso) mso=oso[k];
	}

	// Set up the cumulative index tables
	cu_m=new int[mp+np+op+3];cu_n=cu_m+mp+1;cu_o=cu_n+np+1;
	*cu_m=0;for(l=0;l<mp;l++) cu_m[l+1]=cu_m[l]+osm[l];
	*cu_n=0;for(l=0;l<np;l++) cu_n[l+1]=cu_n[l]+osn[l];
	*cu_o=0;for(l=0;l<op;l++) cu_o[l+1]=cu_o[l]+oso[l];

}

/** Outputs the object boundaries in a binary version of the POV-Ray mesh2
 * format.
 * \param[in] filename the output file to write. */
void fluid_3d::output_contours(const char *filename) {

	// Initialize to suppress warnings
	FILE *outf=NULL;
	std::vector<double> pv;std::vector<int> q,q2;
	MPI_Request sreq[4];

	// Assemble the output filename and open the output file
	if(rank==0) {
		int odata[1];
		*odata=mgmt->n_obj;
		outf=p_safe_fopen(filename,"wb");
		if(fwrite(odata,sizeof(int),1,outf)!=1) {
			fputs("File write error",stderr);
			MPI_Abort(world,1);
		}
	}

	// Output each contour in the list
	for(int i=0; i<mgmt->n_obj; i++) {
		send_contour_data(pv,q,q2,sreq,i+1);
		if(rank==0) output_contour(outf,i);
		MPI_Waitall(4,sreq,MPI_STATUS_IGNORE);
	}

	// Close the output file
	if(rank==0) fclose(outf);
}

/** Outputs a single contour to an open file.
 * \param[in] outf the file handle to write to.
 * \param[in] cnum the number of the contour. */
void fluid_3d::output_contour(FILE *outf,int cnum) {
	int &mp=grid->mp,&np=grid->np,&op=grid->op;
	int l,aip=mp,ajp=np,akp=op,bip=0,bjp=0,bkp=0,cr[4],&ip=*cr,&jp=cr[1],&kp=cr[2];
	int *sz=new int[6*grid->procs], *szp=sz, *szr, *cor=sz + 3*grid->procs, *corp=cor;
	MPI_Request *reqp=reqs;

	// Receive how much data each processors is going to send
	for(kp=0;kp<op;kp++) for(jp=0;jp<np;jp++) for(ip=0;ip<mp;ip++,szp+=3,reqp++) {
		MPI_Cart_rank(grid->cart,cr,cr+3);
		MPI_Irecv(szp,3,MPI_INT,cr[3],fmsg_contour1,grid->cart,reqp);
	}
	MPI_Waitall(grid->procs,reqs,stats);

	// Compute the bounding box of processors involved in this cell, which
	// is used to allocate temporary memory
	int mvert=0,nvert=0,mq2=0,ntri=0,ijk=0;
	for(szp=szr=sz;szp<sz+3*grid->procs;szp+=3,ijk++) {

		// Count global measures of the contour
		if(*szp>mvert) mvert=*szp;
		nvert+=*szp;
		if(szp[1]>mq2) mq2=szp[1];
		ntri+=szp[2];

		// Find bounds of this contour and store a list of processors
		// that have vertices to send
		if(szp[1]>0) {
			*(corp++)=l=ijk%mp;
			if(l<aip) aip=l;
			if(l>bip) bip=l;
			*(corp++)=l=(ijk/mp)%np;
			if(l<ajp) ajp=l;
			if(l>bjp) bjp=l;
			*(corp++)=l=ijk/(mp*np);
			if(l<akp) akp=l;
			if(l>bkp) bkp=l;
			*(szr++)=*szp;
			*(szr++)=szp[1];
			*(szr++)=szp[2];
		}
	}

	int odata[2];*odata=nvert;odata[1]=ntri;
	//printf("# Object %d: nvert=%d ntri=%d\n",cnum,*odata,odata[1]);
	fwrite(odata,sizeof(int),2,outf);

	// Skip completely if there are no vertices
	if(nvert==0) {
		delete [] sz;
		return;
	}

	// Allocate memory for the normal vectors and psoitions
	float *nl=new float[3*nvert],*nlp=nl,*pts=new float[3*mvert],*ptsp;

	// Allocate memory
	double *pinfo=new double[4*mvert];
	int *vt=new int[mvert],*vbx=new int[mq2],*tri=new int[3*ntri],
	    *trip=tri,

	// Compute dimensions of mapping regions, adding one if the last
	// processor in a grid is being used
	    lrm=*max_sizes+1,rm=cu_m[bip+1]-cu_m[aip]+1,
	    lrn=max_sizes[1]+1,rn=cu_n[bjp+1]-cu_n[ajp]+1,
	    lro=max_sizes[2]+1,

	// Allocate memory for mapping regions
	    mrs=lrn*lro+rm*lro+rm*rn,*mrx=new int[4*mrs],
	    *mry=mrx+4*lrn*lro,*mrz=mry+4*rm*lro,*ml1,*ml2,*ml3,*mu1,*mu2,*mu3,
	    *fm=new int[3*lrm*lrn*lro],

	// Begin reading the details of the mesh, storing some of it, and
	// printing the vertex vectors
	    vi,vj,vk,vijk,vn=0,type,*vbxp,*vbxe;
	for(int oo=0;oo<3*lrm*lrn*lro;oo++) fm[oo]=0;
	char *pp;
	double coords[3],&x=*coords,&y=coords[1],&z=coords[2];
	bool posi,posj,posk;
	for(int *szp=sz,*corp=cor;szp<szr;szp+=3,corp+=3) {
		int &i=*corp,&j=corp[1],&k=corp[2],
		    lsm4=osm[i]+4,lsn4=osn[j]+4,lsmn4=lsm4*lsn4;
		MPI_Cart_rank(grid->cart,corp,cr+3);

		ptsp=pts;

		// Set up the boundary layers
		if(i&1) {ml1=mrx;mu1=mrx+lrn*lro;} else {ml1=mrx+lrn*lro;mu1=mrx;}
		if(j&1) {ml2=mry;mu2=mry+rm*lro;} else {ml2=mry+rm*lro;mu2=mry;}
		if(k&1) {ml3=mrz;mu3=mrz+rm*rn;} else {ml3=mrz+rm*rn;mu3=mrz;}

		// Deal with the mesh points if there are any
		if(*szp>0) {

			// Determine if this box borders any of the upper
			// regions of the grid
			posi=i!=mp-1;
			posj=j!=np-1;
			posk=k!=op-1;

			// Receive index, position, and normal information
			MPI_Recv(vt,*szp,MPI_INT,cr[3],fmsg_contour2,grid->cart,MPI_STATUS_IGNORE);
			MPI_Recv(pinfo,4*(*szp),MPI_DOUBLE,cr[3],fmsg_contour3,grid->cart,MPI_STATUS_IGNORE);

			// Loop over the available points
			for(l=0;l<*szp;l++) {

				// Decode the position
				type=vt[l]&3;
				vijk=vt[l]>>2;
				vi=vijk%lsm4;
				vj=(vijk/lsm4)%lsn4;
				vk=vijk/lsmn4;

				// Calculate and store the position
				x=gx[cu_m[i]+vi];
				y=gy[cu_n[j]+vj];
				z=gz[cu_o[k]+vk];
				switch(type) {
					case 0: x+=pinfo[l<<2];break;
					case 1: y+=pinfo[l<<2];break;
					case 2: z+=pinfo[l<<2];break;
				}
				*(ptsp++)=x;*(ptsp++)=y;*(ptsp++)=z;

				// Store the normal vector
				*(nlp++)=pinfo[(l<<2)+1];
				*(nlp++)=pinfo[(l<<2)+2];
				*(nlp++)=pinfo[(l<<2)+3];

				// Store the index for later matching
				if(posk&&vk==oso[k]) {
					if(type==2) p_die();
					mu3[((cu_m[i]-cu_m[aip]+vi+rm*(cu_n[j]-cu_n[ajp]+vj))<<1)+type]=vn;
				}
				if(posj&&vj==osn[j]) {
					if(type==1) p_die();
					mu2[((cu_m[i]-cu_m[aip]+vi+rm*vk)<<1)+(type>>1)]=vn;
				}
				if(posi&&vi==osm[i]) {
					if(type==0) p_die();
					mu1[((vj+lrn*vk)<<1)+type-1]=vn;
				}
				fm[(vi+lrm*(vj+lrn*vk))*3+type]=vn++;
			}

			// Output the positions
			fwrite(pts,sizeof(float),3*(*szp),outf);
		}

		// Reindex the face indices
		if(sz[1]>0) {

			// Determine if this box borders any of the lower
			// regions of the grid
			posi=i==0;
			posj=j==0;
			posk=k==0;

			// Receive the block information
			MPI_Recv(vbx,szp[1],MPI_INT,cr[3],fmsg_contour4,grid->cart,MPI_STATUS_IGNORE);

			vbxp=vbx;vbxe=vbx+szp[1];
			int d1=3*lrm,d2=d1*lrn,key;
			while(vbxp<vbxe) {

				// Get coordinates and then convert into format
				// for looking up edge indices
				vijk=*(vbxp++);
				vi=vijk%lsm4;
				vj=(vijk/lsm4)%lsn4;
				vk=vijk/lsmn4;
				vijk=3*(vi+lrm*(vj+lrn*vk));

				// Look up the indices of the vertices involved
				key=*(vbxp++);
				for(pp=ptri_poly[key];pp<ptri_poly[key+1];pp++,trip++) switch(*pp) {
					case 0: *trip=posk||vk!=0?
						(posj||vj!=0?fm[vijk]:ml2[(cu_m[i]-cu_m[aip]+vi+rm*vk)<<1]):
						ml3[(cu_m[i]-cu_m[aip]+vi+rm*(cu_n[j]-cu_n[ajp]+vj))<<1];break;
					case 1: *trip=posk||vk!=0?fm[vijk+4]:
						ml3[((cu_m[i]-cu_m[aip]+vi+1+rm*(cu_n[j]-cu_n[ajp]+vj))<<1)+1];break;
					case 2: *trip=posk||vk!=0?fm[vijk+d1]:
						ml3[(cu_m[i]-cu_m[aip]+vi+rm*(cu_n[j]-cu_n[ajp]+vj+1))<<1];break;
					case 3: *trip=posk||vk!=0?
						(posi||vi!=0?fm[vijk+1]:ml1[((vj+lrn*vk)<<1)]):
						ml3[((cu_m[i]-cu_m[aip]+vi+rm*(cu_n[j]-cu_n[ajp]+vj))<<1)+1];break;
					case 4: *trip=posj||vj!=0?fm[vijk+d2]:
						ml2[(cu_m[i]-cu_m[aip]+vi+rm*(vk+1))<<1];break;
					case 5: *trip=fm[vijk+d2+4];break;
					case 6: *trip=fm[vijk+d2+d1];break;
					case 7: *trip=posi||vi!=0?fm[vijk+d2+1]:
						ml1[((vj+lrn*(vk+1))<<1)];break;
					case 8: *trip=posj||vj!=0?
						(posi||vi!=0?fm[vijk+2]:ml1[((vj+lrn*vk)<<1)+1]):
						ml2[((cu_m[i]-cu_m[aip]+vi+rm*vk)<<1)+1];break;
					case 9: *trip=posj||vj!=0?fm[vijk+5]:
						ml2[((cu_m[i]-cu_m[aip]+vi+1+rm*vk)<<1)+1];break;
					case 10: *trip=fm[vijk+d1+5];break;
					case 11: *trip=posi||vi!=0?fm[vijk+d1+2]:
						 ml1[((vj+1+lrn*vk)<<1)+1];break;
				}
			}
		}
	}

	// Print the normal vectors
	fwrite(nl,sizeof(float),3*nvert,outf);

	// Print out the face indices
	fwrite(tri,sizeof(int),3*ntri,outf);

	// Deallocate the dynamically allocated memory
	delete [] fm;
	delete [] mrx;
	delete [] tri;
	delete [] vbx;
	delete [] vt;
	delete [] pinfo;
	delete [] pts;
	delete [] nl;
	delete [] sz;
}

void fluid_3d::send_contour_data(std::vector<double> &pv,std::vector<int> &q,std::vector<int> &q2,MPI_Request *sreq,int cnum) {
	pv.clear();q.clear();q2.clear();
	unsigned int ttri=0;
	//geometry &g = *grid;
	bool posi=grid->ip==0,posj=grid->jp==0,posk=grid->kp==0;
	//int ui=(grid->ip==grid->mp-1)?sm-1:sm,uj=(grid->jp==grid->np-1)?sn-1:sn,uk=(grid->kp==grid->op-1)?so-1:so;

	for(int k=0;k<=so;k++) for(int j=0;j<=sn;j++) for(int i=0;i<=sm;i++) {
		int ijk=i+sm4*j+smn4*k;
		if((k>0||posk)&&(j>0||posj)&&i<sm) edge_detect(ijk,0,pv,q,cnum);
		if((k>0||posk)&&j<sn&&(i>0||posi)) edge_detect(ijk,1,pv,q,cnum);
		if(k<so&&(j>0||posj)&&(i>0||posi)) edge_detect(ijk,2,pv,q,cnum);
		if(i<sm&&j<sn&&k<so) box_detect(ijk,q2,ttri,cnum);
		//if(i<ui&&j<uj&&k<uk) box_detect(ijk,q2,ttri,cval);
	}

	// Pack the sizes and send them to the master processor
	*sizes=q.size();
	sizes[1]=q2.size();
	sizes[2]=ttri;
	MPI_Isend(sizes,3,MPI_INT,0,fmsg_contour1,grid->cart,sreq);

	// Send information about vertices to the master processor
	if(sizes[0]>0) {
		MPI_Isend(&q[0],q.size(),MPI_INT,0,fmsg_contour2,grid->cart,sreq+1);
		MPI_Isend(&pv[0],4*q.size(),MPI_DOUBLE,0,fmsg_contour3,grid->cart,sreq+2);
	} else sreq[1]=sreq[2]=MPI_REQUEST_NULL;

	// Send information about blocks to the master processor
	if(sizes[1]>0) MPI_Isend(&q2[0],q2.size(),MPI_INT,0,fmsg_contour4,grid->cart,sreq+3);
	else sreq[3]=MPI_REQUEST_NULL;
}

/* Determines whether the zero level set intersects a given edge, and if so
 * adds it to the buffer to send to the master processor.
 * \param[in] ijk the point to consider.
 * \param[in] dir the direction of the edge (0=x, 1=y, 2=z). */
void fluid_3d::edge_detect(int ijk,int dir,std::vector<double> &pv,std::vector<int> &q,int cnum) {
	const double tol=std::numeric_limits<double>::epsilon();
	int ind=(ijk<<2)|dir;
	int ijk2=ijk+(dir==0?1:(dir==1?sm4:smn4));
	double v1=phi(ijk,cnum),v2=phi(ijk2,cnum);

	if(v1>0?v2<=0:v2>0) {
		double nx,ny,nz,nx2,ny2,nz2;
		double fa=v2-v1,fb=fabs(fa)<tol?0.5:-v1/fa;
		pv.push_back(fb*(dir==0?dx:(dir==1?dy:dz)));
		q.push_back(ind);
		grad_field(ijk,nx,ny,nz,cnum);
		grad_field(ijk2,nx2,ny2,nz2,cnum);
		nx=(1-fb)*nx+fb*nx2;
		ny=(1-fb)*ny+fb*ny2;
		nz=(1-fb)*nz+fb*nz2;
		normalize(nx,ny,nz);
		pv.push_back(nx);
		pv.push_back(ny);
		pv.push_back(nz);
	}
}

/** Computes the gradient of the field at a gridpoint.
 * \param[in] ijk the point to consider.
 * \param[in] (nx,ny,nz) the components of gradient. */
void fluid_3d::grad_field(int ijk,double &nx,double &ny,double &nz,int cnum) {
	nx=dxsp*(phi(ijk+1,cnum)-phi(ijk-1,cnum));
	ny=dysp*(phi(ijk+sm4,cnum)-phi(ijk-sm4,cnum));
	nz=dzsp*(phi(ijk+smn4,cnum)-phi(ijk-smn4,cnum));
}

/* Examines whether a given box is part of a boundary, and if so, adds the
 * index to the list to send to the master processor.
 * \param[in] ijk the index fo the box to consider.
 * \param[in,out] q2 a vector to add the index and key to.
 * \param[in,out] ttri the running total number of triangles. */
void fluid_3d::box_detect(int ijk,std::vector<int> &q2,unsigned int &ttri,int cnum) {

	// Test the eight corners of the box and assemble the key
	int key=phi(ijk,cnum)>0?1:0;
	if(phi(ijk+1,cnum)>0) key|=2;
	if(phi(ijk+sm4+1,cnum)>0) key|=4;
	if(phi(ijk+sm4,cnum)>0) key|=8;
	if(phi(ijk+smn4,cnum)>0) key|=16;
	if(phi(ijk+smn4+1,cnum)>0) key|=32;
	if(phi(ijk+smn4+sm4+1,cnum)>0) key|=64;
	if(phi(ijk+smn4+sm4,cnum)>0) key|=128;

	// If contour passes through this box, then store the key value and box
	// index
	if(key!=255&&key!=0) {
		q2.push_back(ijk);
		q2.push_back(key);
		ttri+=n_poly[key];
	}
}

/** Normalizes a vector.
 * \param[in,out] (nx,ny,nz) the components of the vector, which upon
 *			     completion will have length 1. */
void fluid_3d::normalize(double &nx,double &ny,double &nz) {
	const double tol=std::numeric_limits<double>::epsilon();
	double rsq=nx*nx+ny*ny+nz*nz;
	if(rsq>tol*tol) {
		rsq=1./sqrt(rsq);
		nx*=rsq;
		ny*=rsq;
		nz*=rsq;
	}
}

/** The number of triangles in the contour configuration of each block. */
const char fluid_3d::n_poly[256]={
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,2,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,3,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,3,
	2,3,3,2,3,4,4,3,3,4,4,3,4,5,5,2,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,3,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,4,
	2,3,3,4,3,4,2,3,3,4,4,5,4,5,3,2,
	3,4,4,3,4,5,3,2,4,5,5,4,5,2,4,1,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,3,
	2,3,3,4,3,4,4,5,3,2,4,3,4,3,5,2,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,4,
	3,4,4,3,4,5,5,4,4,3,5,2,5,4,2,1,
	2,3,3,4,3,4,4,5,3,4,4,5,2,3,3,2,
	3,4,4,5,4,5,5,2,4,3,5,4,3,2,4,1,
	3,4,4,5,4,5,3,4,4,5,5,2,3,4,2,1,
	2,3,3,2,3,4,2,1,3,2,4,1,2,1,1,0
};

/** The connectivity of the triangles in the contour configuration of each
 * block. */
const char fluid_3d::tri_poly[2460]={
	0,8,3,0,1,9,1,8,3,9,8,1,1,2,10,0,8,3,1,2,10,9,2,10,0,2,9,2,8,3,
	2,10,8,10,9,8,3,11,2,0,11,2,8,11,0,1,9,0,2,3,11,1,11,2,1,9,11,9,8,11,
	3,10,1,11,10,3,0,10,1,0,8,10,8,11,10,3,9,0,3,11,9,11,10,9,9,8,10,10,8,11,
	4,7,8,4,3,0,7,3,4,0,1,9,8,4,7,4,1,9,4,7,1,7,3,1,1,2,10,8,4,7,
	3,4,7,3,0,4,1,2,10,9,2,10,9,0,2,8,4,7,2,10,9,2,9,7,2,7,3,7,9,4,
	8,4,7,3,11,2,11,4,7,11,2,4,2,0,4,9,0,1,8,4,7,2,3,11,4,7,11,9,4,11,
	9,11,2,9,2,1,3,10,1,3,11,10,7,8,4,1,11,10,1,4,11,1,0,4,7,11,4,4,7,8,
	9,0,11,9,11,10,11,0,3,4,7,11,4,11,9,9,11,10,9,5,4,9,5,4,0,8,3,0,5,4,
	1,5,0,8,5,4,8,3,5,3,1,5,1,2,10,9,5,4,3,0,8,1,2,10,4,9,5,5,2,10,
	5,4,2,4,0,2,2,10,5,3,2,5,3,5,4,3,4,8,9,5,4,2,3,11,0,11,2,0,8,11,
	4,9,5,0,5,4,0,1,5,2,3,11,2,1,5,2,5,8,2,8,11,4,8,5,10,3,11,10,1,3,
	9,5,4,4,9,5,0,8,1,8,10,1,8,11,10,5,4,0,5,0,11,5,11,10,11,0,3,5,4,8,
	5,8,10,10,8,11,9,7,8,5,7,9,9,3,0,9,5,3,5,7,3,0,7,8,0,1,7,1,5,7,
	1,5,3,3,5,7,9,7,8,9,5,7,10,1,2,10,1,2,9,5,0,5,3,0,5,7,3,8,0,2,
	8,2,5,8,5,7,10,5,2,2,10,5,2,5,3,3,5,7,7,9,5,7,8,9,3,11,2,9,5,7,
	9,7,2,9,2,0,2,7,11,2,3,11,0,1,8,1,7,8,1,5,7,11,2,1,11,1,7,7,1,5,
	9,5,8,8,5,7,10,1,3,10,3,11,5,7,0,5,0,9,7,11,0,1,0,10,11,10,0,11,10,0,
	11,0,3,10,5,0,8,0,7,5,7,0,11,10,5,7,11,5,10,6,5,0,8,3,5,10,6,9,0,1,
	5,10,6,1,8,3,1,9,8,5,10,6,1,6,5,2,6,1,1,6,5,1,2,6,3,0,8,9,6,5,
	9,0,6,0,2,6,5,9,8,5,8,2,5,2,6,3,2,8,2,3,11,10,6,5,11,0,8,11,2,0,
	10,6,5,0,1,9,2,3,11,5,10,6,5,10,6,1,9,2,9,11,2,9,8,11,6,3,11,6,5,3,
	5,1,3,0,8,11,0,11,5,0,5,1,5,11,6,3,11,6,0,3,6,0,6,5,0,5,9,6,5,9,
	6,9,11,11,9,8,5,10,6,4,7,8,4,3,0,4,7,3,6,5,10,1,9,0,5,10,6,8,4,7,
	10,6,5,1,9,7,1,7,3,7,9,4,6,1,2,6,5,1,4,7,8,1,2,5,5,2,6,3,0,4,
	3,4,7,8,4,7,9,0,5,0,6,5,0,2,6,7,3,9,7,9,4,3,2,9,5,9,6,2,6,9,
	3,11,2,7,8,4,10,6,5,5,10,6,4,7,2,4,2,0,2,7,11,0,1,9,4,7,8,2,3,11,
	5,10,6,9,2,1,9,11,2,9,4,11,7,11,4,5,10,6,8,4,7,3,11,5,3,5,1,5,11,6,
	5,1,11,5,11,6,1,0,11,7,11,4,0,4,11,0,5,9,0,6,5,0,3,6,11,6,3,8,4,7,
	6,5,9,6,9,11,4,7,9,7,11,9,10,4,9,6,4,10,4,10,6,4,9,10,0,8,3,10,0,1,
	10,6,0,6,4,0,8,3,1,8,1,6,8,6,4,6,1,10,1,4,9,1,2,4,2,6,4,3,0,8,
	1,2,9,2,4,9,2,6,4,0,2,4,4,2,6,8,3,2,8,2,4,4,2,6,10,4,9,10,6,4,
	11,2,3,0,8,2,2,8,11,4,9,10,4,10,6,3,11,2,0,1,6,0,6,4,6,1,10,6,4,1,
	6,1,10,4,8,1,2,1,11,8,11,1,9,6,4,9,3,6,9,1,3,11,6,3,8,11,1,8,1,0,
	11,6,1,9,1,4,6,4,1,3,11,6,3,6,0,0,6,4,6,4,8,11,6,8,7,10,6,7,8,10,
	8,9,10,0,7,3,0,10,7,0,9,10,6,7,10,10,6,7,1,10,7,1,7,8,1,8,0,10,6,7,
	10,7,1,1,7,3,1,2,6,1,6,8,1,8,9,8,6,7,2,6,9,2,9,1,6,7,9,0,9,3,
	7,3,9,7,8,0,7,0,6,6,0,2,7,3,2,6,7,2,2,3,11,10,6,8,10,8,9,8,6,7,
	2,0,7,2,7,11,0,9,7,6,7,10,9,10,7,1,8,0,1,7,8,1,10,7,6,7,10,2,3,11,
	11,2,1,11,1,7,10,6,1,6,7,1,8,9,6,8,6,7,9,1,6,11,6,3,1,3,6,0,9,1,
	11,6,7,7,8,0,7,0,6,3,11,0,11,6,0,7,11,6,7,6,11,3,0,8,11,7,6,0,1,9,
	11,7,6,8,1,9,8,3,1,11,7,6,10,1,2,6,11,7,1,2,10,3,0,8,6,11,7,2,9,0,
	2,10,9,6,11,7,6,11,7,2,10,3,10,8,3,10,9,8,7,2,3,6,2,7,7,0,8,7,6,0,
	6,2,0,2,7,6,2,3,7,0,1,9,1,6,2,1,8,6,1,9,8,8,7,6,10,7,6,10,1,7,
	1,3,7,10,7,6,1,7,10,1,8,7,1,0,8,0,3,7,0,7,10,0,10,9,6,10,7,7,6,10,
	7,10,8,8,10,9,6,8,4,11,8,6,3,6,11,3,0,6,0,4,6,8,6,11,8,4,6,9,0,1,
	9,4,6,9,6,3,9,3,1,11,3,6,6,8,4,6,11,8,2,10,1,1,2,10,3,0,11,0,6,11,
	0,4,6,4,11,8,4,6,11,0,2,9,2,10,9,10,9,3,10,3,2,9,4,3,11,3,6,4,6,3,
	8,2,3,8,4,2,4,6,2,0,4,2,4,6,2,1,9,0,2,3,4,2,4,6,4,3,8,1,9,4,
	1,4,2,2,4,6,8,1,3,8,6,1,8,4,6,6,10,1,10,1,0,10,0,6,6,0,4,4,6,3,
	4,3,8,6,10,3,0,3,9,10,9,3,10,9,4,6,10,4,4,9,5,7,6,11,0,8,3,4,9,5,
	11,7,6,5,0,1,5,4,0,7,6,11,11,7,6,8,3,4,3,5,4,3,1,5,9,5,4,10,1,2,
	7,6,11,6,11,7,1,2,10,0,8,3,4,9,5,7,6,11,5,4,10,4,2,10,4,0,2,3,4,8,
	3,5,4,3,2,5,10,5,2,11,7,6,7,2,3,7,6,2,5,4,9,9,5,4,0,8,6,0,6,2,
	6,8,7,3,6,2,3,7,6,1,5,0,5,4,0,6,2,8,6,8,7,2,1,8,4,8,5,1,5,8,
	9,5,4,10,1,6,1,7,6,1,3,7,1,6,10,1,7,6,1,0,7,8,7,0,9,5,4,4,0,10,
	4,10,5,0,3,10,6,10,7,3,7,10,7,6,10,7,10,8,5,4,10,4,8,10,6,9,5,6,11,9,
	11,8,9,3,6,11,0,6,3,0,5,6,0,9,5,0,11,8,0,5,11,0,1,5,5,6,11,6,11,3,
	6,3,5,5,3,1,1,2,10,9,5,11,9,11,8,11,5,6,0,11,3,0,6,11,0,9,6,5,6,9,
	1,2,10,11,8,5,11,5,6,8,0,5,10,5,2,0,2,5,6,11,3,6,3,5,2,10,3,10,5,3,
	5,8,9,5,2,8,5,6,2,3,8,2,9,5,6,9,6,0,0,6,2,1,5,8,1,8,0,5,6,8,
	3,8,2,6,2,8,1,5,6,2,1,6,1,3,6,1,6,10,3,8,6,5,6,9,8,9,6,10,1,0,
	10,0,6,9,5,0,5,6,0,0,3,8,5,6,10,10,5,6,11,5,10,7,5,11,11,5,10,11,7,5,
	8,3,0,5,11,7,5,10,11,1,9,0,10,7,5,10,11,7,9,8,1,8,3,1,11,1,2,11,7,1,
	7,5,1,0,8,3,1,2,7,1,7,5,7,2,11,9,7,5,9,2,7,9,0,2,2,11,7,7,5,2,
	7,2,11,5,9,2,3,2,8,9,8,2,2,5,10,2,3,5,3,7,5,8,2,0,8,5,2,8,7,5,
	10,2,5,9,0,1,5,10,3,5,3,7,3,10,2,9,8,2,9,2,1,8,7,2,10,2,5,7,5,2,
	1,3,5,3,7,5,0,8,7,0,7,1,1,7,5,9,0,3,9,3,5,5,3,7,9,8,7,5,9,7,
	5,8,4,5,10,8,10,11,8,5,0,4,5,11,0,5,10,11,11,3,0,0,1,9,8,4,10,8,10,11,
	10,4,5,10,11,4,10,4,5,11,3,4,9,4,1,3,1,4,2,5,1,2,8,5,2,11,8,4,5,8,
	0,4,11,0,11,3,4,5,11,2,11,1,5,1,11,0,2,5,0,5,9,2,11,5,4,5,8,11,8,5,
	9,4,5,2,11,3,2,5,10,3,5,2,3,4,5,3,8,4,5,10,2,5,2,4,4,2,0,3,10,2,
	3,5,10,3,8,5,4,5,8,0,1,9,5,10,2,5,2,4,1,9,2,9,4,2,8,4,5,8,5,3,
	3,5,1,0,4,5,1,0,5,8,4,5,8,5,3,9,0,5,0,3,5,9,4,5,4,11,7,4,9,11,
	9,10,11,0,8,3,4,9,7,9,11,7,9,10,11,1,10,11,1,11,4,1,4,0,7,4,11,3,1,4,
	3,4,8,1,10,4,7,4,11,10,11,4,4,11,7,9,11,4,9,2,11,9,1,2,9,7,4,9,11,7,
	9,1,11,2,11,1,0,8,3,11,7,4,11,4,2,2,4,0,11,7,4,11,4,2,8,3,4,3,2,4,
	2,9,10,2,7,9,2,3,7,7,4,9,9,10,7,9,7,4,10,2,7,8,7,0,2,0,7,3,7,10,
	3,10,2,7,4,10,1,10,0,4,0,10,1,10,2,8,7,4,4,9,1,4,1,7,7,1,3,4,9,1,
	4,1,7,0,8,1,8,7,1,4,0,3,7,4,3,4,8,7,9,10,8,10,11,8,3,0,9,3,9,11,
	11,9,10,0,1,10,0,10,8,8,10,11,3,1,10,11,3,10,1,2,11,1,11,9,9,11,8,3,0,9,
	3,9,11,1,2,9,2,11,9,0,2,11,8,0,11,3,2,11,2,3,8,2,8,10,10,8,9,9,10,2,
	0,9,2,2,3,8,2,8,10,0,1,8,1,10,8,1,10,2,1,3,8,9,1,8,0,9,1,0,3,8
};


/** Writes per-frame 3D fields used for budget and spectral energy-exchange post-processing.
 *
 * The files are native-endian binary files with this layout:
 *   8 bytes:  ASCII magic string "RMT3D3D\0"
 *   3 int32:  global dimensions m, n, o
 *   m*n*o float32 values, with x-index varying fastest.
 *
 * The legacy ppe3d/ppf3d definitions are intentionally replaced here:
 *   ppe3d = 2*nu*chi_f*S'_{ij}S'_{ij}
 *   ppf3d = chi_f*u'_fluid_i*fp_total_i, where fp_total is the particle
 *           feedback acceleration cached from the actual momentum RHS during
 *           the velocity update; ppf3d_full is the full-domain counterpart.
 *
 * HDF5-style dataset definitions and attributes for these native binary files
 * are written to frame_fields.<frame>.meta so downstream converters can create
 * the requested /fields, /diagnostics, and /particles HDF5 hierarchy without
 * ambiguity.
 */
namespace {
    double g_frame_ubar_full[3] = {0.,0.,0.};
    double g_frame_ubar_fluid[3] = {0.,0.,0.};
    const char *g_deriv_op = "solver-consistent face-normal velocity gradients averaged to cell centers (same face-difference convention as fluid_stress); fp_total is cached from momentum RHS";

    void write_attr(FILE *fh,const char *path,const char *key,const char *value) {
        fprintf(fh, "attr %s %s = %s\n", path, key, value);
    }

    void write_scalar_attr(FILE *fh,const char *path,const char *definition,const char *units,
                           const char *sign,const fluid_3d *f3d,int frame_num) {
        write_attr(fh,path,"units",units);
        write_attr(fh,path,"definition",definition);
        write_attr(fh,path,"sign_convention",sign);
        write_attr(fh,path,"normalization","pointwise value on the solver cell-centered grid; global means are arithmetic spatial means");
        write_attr(fh,path,"mask_used","see definition");
        write_attr(fh,path,"derivative_operator",g_deriv_op);
        fprintf(fh,"attr %s density_used = %.17g\n", path, f3d->mgmt->fm.rho);
        fprintf(fh,"attr %s viscosity_used = %.17g\n", path, f3d->mgmt->fm.mu/f3d->mgmt->fm.rho);
        fprintf(fh,"attr %s time = %.17g\n", path, f3d->time);
        fprintf(fh,"attr %s timestep = %d\n", path, f3d->nt);
        fprintf(fh,"attr %s frame_index = %d\n", path, frame_num);
    }
}

void fluid_3d::compute_frame_ubar(double (&ubar_full)[3],double (&ubar_fluid)[3]) {
    double local_full[3] = {0.,0.,0.};
    double local_fluid[4] = {0.,0.,0.,0.};
    for(int k=0;k<so;k++) for(int j=0;j<sn;j++) for(int i=0;i<sm;i++) {
        const int eid = index(i,j,k);
        field &f = u0[eid];
        const double chi_f = chi_f_cell_value(eid);
        for(int c=0;c<3;c++) {
            local_full[c] += f.vel[c];
            local_fluid[c] += chi_f*f.vel[c];
        }
        local_fluid[3] += chi_f;
    }
    double global_full[3] = {0.,0.,0.};
    double global_fluid[4] = {0.,0.,0.,0.};
    MPI_Allreduce(local_full, global_full, 3, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(local_fluid, global_fluid, 4, MPI_DOUBLE, MPI_SUM, world);
    const double denom = static_cast<double>(m)*static_cast<double>(n)*static_cast<double>(o);
    for(int c=0;c<3;c++) {
        ubar_full[c] = global_full[c]/denom;
        ubar_fluid[c] = global_fluid[c]/std::max(global_fluid[3], 1.e-30);
    }
}

void fluid_3d::write_3d_exchange_fields(const char* dirname, const int frame_num) {
    compute_frame_ubar(g_frame_ubar_full, g_frame_ubar_fluid);
    compute_stress(false);

    const char* names[] = {
        "u_x", "u_y", "u_z",
        "vx3d", "vy3d", "vz3d",
        "chi_f", "chi_p",
        "ppe3d", "eps_all_3d", "ppf3d", "ppf3d_full",
        "fp_total_x", "fp_total_y", "fp_total_z",
        "fp_stressdiv_x", "fp_stressdiv_y", "fp_stressdiv_z",
        "sigma_total_xx", "sigma_total_xy", "sigma_total_xz", "sigma_total_yy", "sigma_total_yz", "sigma_total_zz",
        "ppe_centered_debug"
    };
    const int field_ids[] = {
        0, 1, 2,
        0, 1, 2,
        3, 4,
        5, 6, 7, 8,
        9, 10, 11,
        12, 13, 14,
        15, 16, 17, 18, 19, 20,
        21
    };
    const int nfields = sizeof(names)/sizeof(names[0]);
    char filename[512];

    for(int f=0; f<nfields; f++) {
        sprintf(filename, "%s/%s.%05d", dirname, names[f], frame_num);
        write_scalar_3d_binary(filename, field_ids[f]);
    }
    write_frame_diagnostics(dirname, frame_num);
    write_3d_metadata(dirname, frame_num);
}

void fluid_3d::write_scalar_3d_binary(const char* filename, int field_type) {
    const int local_len = sm*sn*so;
    std::vector<float> local(local_len);

    int lp = 0;
    for(int k=0; k<so; k++) for(int j=0; j<sn; j++) for(int i=0; i<sm; i++, lp++) {
        const int eid = index(i,j,k);
        local[lp] = static_cast<float>(output_field_3d_value(field_type, eid));
    }

    const int TAG_INFO = 6001;
    const int TAG_DATA = 6002;
    if(rank==0) {
        const size_t global_len = static_cast<size_t>(m)*static_cast<size_t>(n)*static_cast<size_t>(o);
        std::vector<float> global(global_len, 0.f);

        int info[6] = {ai, aj, ak, sm, sn, so};
        for(int kk=0; kk<info[5]; kk++) for(int jj=0; jj<info[4]; jj++) for(int ii=0; ii<info[3]; ii++) {
            const size_t lidx = static_cast<size_t>(ii) + static_cast<size_t>(info[3])*(static_cast<size_t>(jj) + static_cast<size_t>(info[4])*static_cast<size_t>(kk));
            const size_t gidx = static_cast<size_t>(info[0]+ii) + static_cast<size_t>(m)*(static_cast<size_t>(info[1]+jj) + static_cast<size_t>(n)*static_cast<size_t>(info[2]+kk));
            global[gidx] = local[lidx];
        }

        for(int sender=1; sender<grid->procs; sender++) {
            MPI_Recv(info, 6, MPI_INT, sender, TAG_INFO, grid->cart, MPI_STATUS_IGNORE);
            const int recv_len = info[3]*info[4]*info[5];
            std::vector<float> recv_buf(recv_len);
            if(recv_len>0) MPI_Recv(&recv_buf[0], recv_len, MPI_FLOAT, sender, TAG_DATA, grid->cart, MPI_STATUS_IGNORE);

            for(int kk=0; kk<info[5]; kk++) for(int jj=0; jj<info[4]; jj++) for(int ii=0; ii<info[3]; ii++) {
                const size_t lidx = static_cast<size_t>(ii) + static_cast<size_t>(info[3])*(static_cast<size_t>(jj) + static_cast<size_t>(info[4])*static_cast<size_t>(kk));
                const size_t gidx = static_cast<size_t>(info[0]+ii) + static_cast<size_t>(m)*(static_cast<size_t>(info[1]+jj) + static_cast<size_t>(n)*static_cast<size_t>(info[2]+kk));
                global[gidx] = recv_buf[lidx];
            }
        }

        FILE *fh = p_safe_fopen(filename, "wb");
        const char magic[8] = {'R','M','T','3','D','3','D','\0'};
        const int dims[3] = {m, n, o};
        fwrite(magic, sizeof(char), 8, fh);
        fwrite(dims, sizeof(int), 3, fh);
        if(global_len>0) fwrite(&global[0], sizeof(float), global_len, fh);
        fclose(fh);
    } else {
        int info[6] = {ai, aj, ak, sm, sn, so};
        MPI_Send(info, 6, MPI_INT, 0, TAG_INFO, grid->cart);
        if(local_len>0) MPI_Send(&local[0], local_len, MPI_FLOAT, 0, TAG_DATA, grid->cart);
    }

    MPI_Barrier(grid->cart);
}

double fluid_3d::output_field_3d_value(int field_type,int eid) {
    double fp[3] = {0.,0.,0.};
    switch(field_type) {
        case 0: return u0[eid].vel[0];
        case 1: return u0[eid].vel[1];
        case 2: return u0[eid].vel[2];
        case 3: return chi_f_cell_value(eid);
        case 4: return chi_p_cell_value(eid);
        case 5: return ppe_cell_value(eid);
        case 6: return 2.*(mgmt->fm.mu/mgmt->fm.rho)*strain_s2_cell_value(eid);
        case 7: return ppf_cell_value(eid);
        case 8: return ppf_full_cell_value(eid);
        case 9: fp_total_cell_value(eid, fp); return fp[0];
        case 10: fp_total_cell_value(eid, fp); return fp[1];
        case 11: fp_total_cell_value(eid, fp); return fp[2];
        case 12: fp_stressdiv_cell_value(eid, fp); return fp[0];
        case 13: fp_stressdiv_cell_value(eid, fp); return fp[1];
        case 14: fp_stressdiv_cell_value(eid, fp); return fp[2];
        case 15: return 0.5*(u0[eid].sigma[0][0] + u0[eid+1].sigma[0][0]);
        case 16: return 0.25*(u0[eid].sigma[0][1] + u0[eid+1].sigma[0][1] + u0[eid].sigma[1][0] + u0[eid+sm4].sigma[1][0]);
        case 17: return 0.25*(u0[eid].sigma[0][2] + u0[eid+1].sigma[0][2] + u0[eid].sigma[2][0] + u0[eid+smn4].sigma[2][0]);
        case 18: return 0.5*(u0[eid].sigma[1][1] + u0[eid+sm4].sigma[1][1]);
        case 19: return 0.25*(u0[eid].sigma[1][2] + u0[eid+sm4].sigma[1][2] + u0[eid].sigma[2][1] + u0[eid+smn4].sigma[2][1]);
        case 20: return 0.5*(u0[eid].sigma[2][2] + u0[eid+smn4].sigma[2][2]);
        case 21: return 2.*(mgmt->fm.mu/mgmt->fm.rho)*chi_f_cell_value(eid)*centered_strain_s2_cell_value(eid);
        default: return 0.;
    }
}

double fluid_3d::chi_p_cell_value(int eid) {
    return (min_phi(eid) <= 0.) ? 1. : 0.;
}

double fluid_3d::chi_f_cell_value(int eid) {
    return 1. - chi_p_cell_value(eid);
}

double fluid_3d::strain_s2_cell_value(int eid) {
    const int strides[3] = {1, sm4, smn4};
    const double invh[3] = {dxsp, dysp, dzsp};
    double grad[3][3];

    // Solver-consistent gradient: reuse the face-normal difference convention
    // used by fluid_stress(), then average lower/upper face gradients back to
    // the cell center before forming the symmetric strain-rate tensor.
    for(int r=0; r<3; r++) {
        const int step = strides[r];
        for(int c=0; c<3; c++) {
            const double lower = invh[r]*(u0[eid].vel[c] - u0[eid-step].vel[c]);
            const double upper = invh[r]*(u0[eid+step].vel[c] - u0[eid].vel[c]);
            grad[r][c] = 0.5*(lower + upper);
        }
    }

    double s2 = 0.;
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
        const double sij = 0.5*(grad[i][j] + grad[j][i]);
        s2 += sij*sij;
    }
    return s2;
}

double fluid_3d::centered_strain_s2_cell_value(int eid) {
    const int strides[3] = {1, sm4, smn4};
    const double invh[3] = {dxsp, dysp, dzsp};
    double grad[3][3];
    for(int r=0; r<3; r++) {
        const int step = strides[r];
        for(int c=0; c<3; c++) grad[r][c] = 0.5*invh[r]*(u0[eid+step].vel[c] - u0[eid-step].vel[c]);
    }
    double s2 = 0.;
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
        const double sij = 0.5*(grad[i][j] + grad[j][i]);
        s2 += sij*sij;
    }
    return s2;
}

void fluid_3d::fp_stressdiv_cell_value(int eid,double (&fp)[3]) {
    const int strides[3] = {1, sm4, smn4};
    const double invh[3] = {dxsp, dysp, dzsp};
    for(int c=0;c<3;c++) fp[c] = 0.;

    for(int face=0; face<3; face++) {
        double lower_fluid[3] = {0.,0.,0.};
        double upper_fluid[3] = {0.,0.,0.};
        const int strd = strides[face];
        for(int c=0;c<3;c++) {
            lower_fluid[c] = mgmt->fm.mu * invh[face] * (u0[eid].vel[c] - u0[eid-strd].vel[c]);
            upper_fluid[c] = mgmt->fm.mu * invh[face] * (u0[eid+strd].vel[c] - u0[eid].vel[c]);
            const double lower_particle = u0[eid].sigma[face][c] - lower_fluid[c];
            const double upper_particle = u0[eid+strd].sigma[face][c] - upper_fluid[c];
            fp[c] += invh[face]*(upper_particle - lower_particle);
        }
    }
#if defined(VAR_DEN)
    const double rho = lrho[eid+G0];
    if(rho!=0.) for(int c=0;c<3;c++) fp[c] /= rho;
#endif
}

void fluid_3d::fp_total_cell_value(int eid,double (&fp)[3]) {
    for(int c=0;c<3;c++) fp[c] = fp_rhs_total0[3*eid+c];
}

double fluid_3d::ppe_cell_value(int eid) {
    return 2.*(mgmt->fm.mu/mgmt->fm.rho)*chi_f_cell_value(eid)*strain_s2_cell_value(eid);
}

double fluid_3d::ppf_cell_value(int eid) {
    double fp[3] = {0.,0.,0.};
    fp_total_cell_value(eid, fp);
    double val = 0.;
    for(int c=0;c<3;c++) val += (u0[eid].vel[c] - g_frame_ubar_fluid[c])*fp[c];
    return chi_f_cell_value(eid)*val;
}

double fluid_3d::ppf_full_cell_value(int eid) {
    double fp[3] = {0.,0.,0.};
    fp_total_cell_value(eid, fp);
    double val = 0.;
    for(int c=0;c<3;c++) val += (u0[eid].vel[c] - g_frame_ubar_full[c])*fp[c];
    return val;
}

void fluid_3d::write_frame_diagnostics(const char* dirname,const int frame_num) {
    const int pobj = mgmt->n_obj;
    const int pnq = 12; // ppf total/abs/pos/neg, near eps/S2/ppf/abs/pos/neg/count/vol
    const double dvol = dx*dy*dz;
    const double small = 1.e-30;
    const double alpha_near = 1.5, gap_g0 = 1.e-6, gap_pow = 1.;

    // Build true connected-component clusters from near-neighbor graph. This is
    // replicated on every rank because object centers are globally available.
    std::vector<int> near_count(pobj, 0), parent(pobj, 0), cluster_id(pobj, 0), cluster_size(pobj, 1), state(pobj, 0);
    std::vector<double> gap_weight_sum(pobj, 0.);
    for(int i=0;i<pobj;i++) parent[i] = i;
    auto find_root = [&](int a) -> int { while(parent[a]!=a){ parent[a]=parent[parent[a]]; a=parent[a]; } return a; };
    auto unite = [&](int a,int b){ int ra=find_root(a), rb=find_root(b); if(ra!=rb) parent[rb]=ra; };
    double I_gap_total = 0., near_pair_weight_sum = 0.;
    int near_pair_count = 0;
    for(int i=0;i<pobj;i++) for(int j=i+1;j<pobj;j++) {
        object *oi = mgmt->objs[i], *oj = mgmt->objs[j];
        double ddx = oi->c[0] - oj->c[0], ddy = oi->c[1] - oj->c[1], ddz = oi->c[2] - oj->c[2];
        if(mgmt->x_prd) ddx -= mgmt->lx*floor(ddx/mgmt->lx + 0.5);
        if(mgmt->y_prd) ddy -= mgmt->ly*floor(ddy/mgmt->ly + 0.5);
        if(mgmt->z_prd) ddz -= mgmt->lz*floor(ddz/mgmt->lz + 0.5);
        const double rij = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
        const double asum = oi->primary_dim + oj->primary_dim;
        const double dp = 2.*std::max(oi->primary_dim, 1.e-30);
        const double gij = std::max(rij - asum, 0.)/dp;
        const double wij = 1./pow(gij + gap_g0, gap_pow);
        gap_weight_sum[i] += wij;
        gap_weight_sum[j] += wij;
        I_gap_total += wij;
        if(rij < alpha_near*asum) {
            near_count[i]++;
            near_count[j]++;
            near_pair_count++;
            near_pair_weight_sum += wij;
            unite(i,j);
        }
    }
    std::vector<int> comp_size(pobj, 0);
    for(int i=0;i<pobj;i++) comp_size[find_root(i)]++;
    int cluster_count = 0, largest_cluster_size = 0, clustered_particles = 0;
    double mean_cluster_size = 0.;
    for(int i=0;i<pobj;i++) if(parent[i]==i) {
        cluster_count++;
        largest_cluster_size = std::max(largest_cluster_size, comp_size[i]);
        mean_cluster_size += comp_size[i];
    }
    mean_cluster_size = cluster_count>0 ? mean_cluster_size/cluster_count : 0.;
    for(int i=0;i<pobj;i++) {
        const int root = find_root(i);
        cluster_id[i] = root + 1;
        cluster_size[i] = comp_size[root];
        state[i] = (cluster_size[i]==1) ? 0 : ((cluster_size[i]==2) ? 1 : 2);
        if(cluster_size[i]>=3) clustered_particles++;
    }

    double local[40] = {0.};
    local[20] = 1.e300; // min ppe
    std::vector<double> plocal(pobj*pnq, 0.), pglobal(pobj*pnq, 0.);

    for(int k=0;k<so;k++) for(int j=0;j<sn;j++) for(int i=0;i<sm;i++) {
        const int eid = index(i,j,k);
        const double chi_f = chi_f_cell_value(eid), chi_p = 1.-chi_f;
        const double up_full[3] = {u0[eid].vel[0] - g_frame_ubar_full[0], u0[eid].vel[1] - g_frame_ubar_full[1], u0[eid].vel[2] - g_frame_ubar_full[2]};
        const double up_fluid[3] = {u0[eid].vel[0] - g_frame_ubar_fluid[0], u0[eid].vel[1] - g_frame_ubar_fluid[1], u0[eid].vel[2] - g_frame_ubar_fluid[2]};
        const double s2 = strain_s2_cell_value(eid);
        const double eps = ppe_cell_value(eid);
        const double eps_all = 2.*(mgmt->fm.mu/mgmt->fm.rho)*s2;
        const double ppe_centered = 2.*(mgmt->fm.mu/mgmt->fm.rho)*chi_f*centered_strain_s2_cell_value(eid);
        const double ppf_fluid = ppf_cell_value(eid);
        const double ppf_full = ppf_full_cell_value(eid);
        double fp_rhs[3], fp_sd[3]; fp_total_cell_value(eid, fp_rhs); fp_stressdiv_cell_value(eid, fp_sd);
        double fp_rhs_mag2 = 0., fp_diff_mag2 = 0., fp_rhs_linf = 0., fp_diff_linf = 0.;
        for(int c=0;c<3;c++) {
            const double diff = fp_rhs[c]-fp_sd[c];
            fp_rhs_mag2 += fp_rhs[c]*fp_rhs[c];
            fp_diff_mag2 += diff*diff;
            fp_rhs_linf = std::max(fp_rhs_linf, fabs(fp_rhs[c]));
            fp_diff_linf = std::max(fp_diff_linf, fabs(diff));
        }

        local[0] += 0.5*(up_full[0]*up_full[0] + up_full[1]*up_full[1] + up_full[2]*up_full[2]);
        local[1] += 0.5*chi_f*(up_fluid[0]*up_fluid[0] + up_fluid[1]*up_fluid[1] + up_fluid[2]*up_fluid[2]);
        local[2] += eps_all;
        local[3] += eps;
        local[4] += ppf_full;
        local[5] += ppf_full>0. ? ppf_full : 0.;
        local[6] += ppf_full<0. ? -ppf_full : 0.;
        local[7] += fabs(ppf_full);
        local[8] += ppf_fluid;
        local[9] += ppf_fluid>0. ? ppf_fluid : 0.;
        local[10] += ppf_fluid<0. ? -ppf_fluid : 0.;
        local[11] += fabs(ppf_fluid);
        local[12] += fp_rhs_mag2;
        local[13] += fp_diff_mag2;
        local[14] = std::max(local[14], fp_rhs_linf);
        local[15] = std::max(local[15], fp_diff_linf);
        local[16] += (ppe_centered - eps)*(ppe_centered - eps);
        local[17] += eps*eps;
        local[18] += fabs(chi_f + chi_p - 1.);
        local[19] += chi_p;
        if(eps < local[20]) local[20] = eps;

        int assign = -1;
        if(pobj>0 && (fp_rhs_mag2>small || chi_p>0.5)) {
            double best = 1.e300;
            for(int po=0; po<pobj; po++) {
                object *obj = mgmt->objs[po];
                double ddx = lx0[i] - obj->c[0], ddy = ly0[j] - obj->c[1], ddz = lz0[k] - obj->c[2];
                if(mgmt->x_prd) ddx -= mgmt->lx*floor(ddx/mgmt->lx + 0.5);
                if(mgmt->y_prd) ddy -= mgmt->ly*floor(ddy/mgmt->ly + 0.5);
                if(mgmt->z_prd) ddz -= mgmt->lz*floor(ddz/mgmt->lz + 0.5);
                const double r2 = ddx*ddx + ddy*ddy + ddz*ddz;
                if(r2<best) {best=r2; assign=po;}
            }
        }
        if(assign>=0) {
            double *q = &plocal[assign*pnq];
            q[0] += ppf_full*dvol;
            q[1] += fabs(ppf_full)*dvol;
            q[2] += (ppf_full>0. ? ppf_full : 0.)*dvol;
            q[3] += (ppf_full<0. ? -ppf_full : 0.)*dvol;
        }
        for(int po=0; po<pobj; po++) {
            object *obj = mgmt->objs[po];
            double ddx = lx0[i] - obj->c[0], ddy = ly0[j] - obj->c[1], ddz = lz0[k] - obj->c[2];
            if(mgmt->x_prd) ddx -= mgmt->lx*floor(ddx/mgmt->lx + 0.5);
            if(mgmt->y_prd) ddy -= mgmt->ly*floor(ddy/mgmt->ly + 0.5);
            if(mgmt->z_prd) ddz -= mgmt->lz*floor(ddz/mgmt->lz + 0.5);
            const double rnear = 3.0*obj->primary_dim; // 1.5*dp for spherical primary_dim=radius
            if(ddx*ddx + ddy*ddy + ddz*ddz < rnear*rnear) {
                double *q = &plocal[po*pnq];
                q[4] += eps*dvol;
                q[5] += s2*dvol;
                q[6] += ppf_fluid*dvol;
                q[7] += fabs(ppf_fluid)*dvol;
                q[8] += (ppf_fluid>0. ? ppf_fluid : 0.)*dvol;
                q[9] += (ppf_fluid<0. ? -ppf_fluid : 0.)*dvol;
                q[10] += 1.;
                q[11] += dvol;
            }
        }
    }
    if(pobj>0) MPI_Reduce(&plocal[0], &pglobal[0], pobj*pnq, MPI_DOUBLE, MPI_SUM, 0, world);

    double global[40] = {0.};
    MPI_Allreduce(local, global, 14, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(local+14, global+14, 2, MPI_DOUBLE, MPI_MAX, world);
    MPI_Allreduce(local+16, global+16, 4, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(local+20, global+20, 1, MPI_DOUBLE, MPI_MIN, world);

    const double denom = static_cast<double>(m)*static_cast<double>(n)*static_cast<double>(o);
    for(int q=0;q<=13;q++) global[q] /= denom;
    for(int q=16;q<=19;q++) global[q] /= denom;
    const double PPF_full_asym = (global[5]-global[6]) / std::max(global[5]+global[6], small);
    const double PPF_fluid_asym = (global[9]-global[10]) / std::max(global[9]+global[10], small);
    const double fp_l2_rel = sqrt(global[13]) / std::max(sqrt(global[12]), small);
    const double fp_linf_rel = global[15] / std::max(global[14], small);
    const double ppe_centered_rel = sqrt(global[16]) / std::max(sqrt(global[17]), small);
    const double nanv = std::numeric_limits<double>::quiet_NaN();

    if(rank==0) {
        char filename[512];
        sprintf(filename, "%s/frame_diagnostics.%05d", dirname, frame_num);
        FILE *fh = p_safe_fopen(filename, "w");
        fprintf(fh,"K_full %.17g\nK_fluid %.17g\n", global[0], global[1]);
        fprintf(fh,"eps_full %.17g\neps_fluid %.17g\n", global[2], global[3]);
        fprintf(fh,"Wp_full %.17g\nWp_full_pos %.17g\nWp_full_neg %.17g\nWp_full_abs %.17g\nPPF_full_asymmetry %.17g\n", global[4], global[5], global[6], global[7], PPF_full_asym);
        fprintf(fh,"Wp_fluid %.17g\nWp_fluid_pos %.17g\nWp_fluid_neg %.17g\nWp_fluid_abs %.17g\nPPF_fluid_asymmetry %.17g\n", global[8], global[9], global[10], global[11], PPF_fluid_asym);
        fprintf(fh,"ubar_full_x %.17g\nubar_full_y %.17g\nubar_full_z %.17g\n", g_frame_ubar_full[0], g_frame_ubar_full[1], g_frame_ubar_full[2]);
        fprintf(fh,"ubar_fluid_x %.17g\nubar_fluid_y %.17g\nubar_fluid_z %.17g\n", g_frame_ubar_fluid[0], g_frame_ubar_fluid[1], g_frame_ubar_fluid[2]);
        fprintf(fh,"dK_full_dt %.17g\ndK_fluid_dt %.17g\n", nanv, nanv);
        fprintf(fh,"residual_full %.17g\nresidual_full_normalized %.17g\nresidual_fluid %.17g\nresidual_fluid_normalized %.17g\n", nanv, nanv, nanv, nanv);
        fprintf(fh,"fp_rhs_stressdiv_l2_rel_error %.17g\nfp_rhs_stressdiv_linf_rel_error %.17g\n", fp_l2_rel, fp_linf_rel);
        fprintf(fh,"ppe_centered_vs_solver_l2_rel_error %.17g\n", ppe_centered_rel);
        fprintf(fh,"ppf_force_mean %.17g\nppf_stress_mean %.17g\nppf_force_stress_rel_error %.17g\n", global[8], nanv, nanv);
        fprintf(fh,"min_ppe3d %.17g\nchi_sum_mismatch_mean %.17g\nparticle_volume_fraction_cells %.17g\n", global[20], global[18], global[19]);
        fprintf(fh,"cluster_count %d\nlargest_cluster_size %d\nlargest_cluster_fraction %.17g\nmean_cluster_size %.17g\nclustered_particle_fraction %.17g\n", cluster_count, largest_cluster_size, pobj>0?double(largest_cluster_size)/pobj:0., mean_cluster_size, pobj>0?double(clustered_particles)/pobj:0.);
        fprintf(fh,"I_gap_total %.17g\nI_gap_per_particle %.17g\nmean_pair_weight %.17g\n", I_gap_total, pobj>0?I_gap_total/pobj:0., near_pair_count>0?near_pair_weight_sum/near_pair_count:0.);

        double ppf_abs_sum = 0.;
        std::vector<double> ppf_abs(pobj, 0.);
        for(int p=0;p<pobj;p++) { ppf_abs[p] = pglobal[p*pnq+1]; ppf_abs_sum += ppf_abs[p]; }
        std::vector<double> sorted_abs = ppf_abs;
        std::sort(sorted_abs.begin(), sorted_abs.end(), [](double a,double b){return a>b;});
        auto top_share = [&](double frac) -> double { int nsel = pobj>0 ? std::max(1, int(ceil(frac*pobj))) : 0; double s=0.; for(int i=0;i<nsel && i<pobj;i++) s+=sorted_abs[i]; return s/std::max(ppf_abs_sum, small); };
        double hhi=0., cluster_share=0., pair_share=0., iso_share=0.;
        for(int p=0;p<pobj;p++) {
            const double w = ppf_abs[p]/std::max(ppf_abs_sum, small);
            hhi += w*w;
            if(state[p]==2) cluster_share += ppf_abs[p];
            else if(state[p]==1) pair_share += ppf_abs[p];
            else iso_share += ppf_abs[p];
        }
        fprintf(fh,"PPF_top10_share %.17g\nPPF_top20_share %.17g\nPPF_Herfindahl_index %.17g\nPPF_cluster_share %.17g\nPPF_pair_like_share %.17g\nPPF_isolated_share %.17g\n", top_share(0.10), top_share(0.20), hhi, cluster_share/std::max(ppf_abs_sum, small), pair_share/std::max(ppf_abs_sum, small), iso_share/std::max(ppf_abs_sum, small));
        fclose(fh);

        sprintf(filename, "%s/particles.%05d", dirname, frame_num);
        fh = p_safe_fopen(filename, "w");
        fprintf(fh,"# id center_x center_y center_z radius equivalent_radius volume deformation_available aspect_ratio principal_axis_x principal_axis_y principal_axis_z near_neighbor_count gap_weight_sum cluster_id cluster_size state ppf_i_total ppf_i_abs ppf_i_pos ppf_i_neg eps_near_i S2_near_i ppf_near_i_total ppf_near_i_abs ppf_near_i_pos ppf_near_i_neg local_neighborhood_cell_count local_neighborhood_volume\n");
        for(int p=0;p<pobj;p++) {
            object *obj = mgmt->objs[p];
            const double *q = &pglobal[p*pnq];
            const double vnear = std::max(q[11], small);
            fprintf(fh,"%d %.17g %.17g %.17g %.17g %.17g %.17g %d %.17g %.17g %.17g %.17g %d %.17g %d %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                    p+1, obj->c[0], obj->c[1], obj->c[2], obj->primary_dim, obj->primary_dim, obj->volume,
                    0, nanv, nanv, nanv, nanv, near_count[p], gap_weight_sum[p], cluster_id[p], cluster_size[p], state[p],
                    q[0], q[1], q[2], q[3], q[4]/vnear, q[5]/vnear, q[6]/vnear, q[7]/vnear, q[8]/vnear, q[9]/vnear, q[10], q[11]);
        }
        fclose(fh);
    }
}

void fluid_3d::write_3d_metadata(const char* dirname,const int frame_num) {
    if(rank!=0) return;
    char filename[512];
    sprintf(filename, "%s/frame_fields.%05d.meta", dirname, frame_num);
    FILE *fh = p_safe_fopen(filename, "w");
    fprintf(fh,"# Native binary datasets mirror the requested HDF5 paths. Convert name.<frame> to /fields/name.\n");
    fprintf(fh,"dimensions %d %d %d\n", m,n,o);
    fprintf(fh,"diagnostics_file frame_diagnostics.%05d\n", frame_num);
    fprintf(fh,"particles_file particles.%05d\n", frame_num);
    write_scalar_attr(fh,"/fields/ppe3d","2 * nu * chi_f * S'_{ij}S'_{ij}, where S' is computed from velocity fluctuations using the solver-consistent face-normal fluid_stress operator averaged to cell centers","length^2/time^3","positive definite carrier-fluid viscous dissipation density per unit mass",this,frame_num);
    write_scalar_attr(fh,"/fields/eps_all_3d","2 * nu * S'_{ij}S'_{ij}; full-domain one-fluid diagnostic using the solver-consistent face-normal fluid_stress operator averaged to cell centers","length^2/time^3","non-negative diagnostic",this,frame_num);
    write_scalar_attr(fh,"/fields/ppe_centered_debug","debug-only centered-difference 2 * nu * chi_f * S_{ij}S_{ij}; not the formal budget ppe3d unless identical to the solver operator","length^2/time^3","debug non-negative diagnostic",this,frame_num);
    write_scalar_attr(fh,"/fields/ppf3d","carrier-fluid masked power density chi_f * u'_fluid_i * fp_total_i, where fp_total_i is cached from the actual particle-feedback acceleration used in the momentum RHS","length^2/time^3","positive means particle feedback injects kinetic energy into carrier-fluid velocity fluctuation; negative means extraction",this,frame_num);
    write_scalar_attr(fh,"/fields/ppf3d_full","full-domain one-fluid power density u'_full_i * fp_total_i for spectral velocity-force cross budgets","length^2/time^3","positive means particle feedback injects kinetic energy into full-domain velocity fluctuation; negative means extraction",this,frame_num);
    write_scalar_attr(fh,"/fields/fp_total_x","fp_total is the actual particle-feedback acceleration used in the momentum RHS; x component cached during compute_ustar/acceleration, not reconstructed at output time","acceleration","positive x acceleration on the momentum RHS",this,frame_num);
    write_scalar_attr(fh,"/fields/fp_total_y","fp_total is the actual particle-feedback acceleration used in the momentum RHS; y component cached during compute_ustar/acceleration, not reconstructed at output time","acceleration","positive y acceleration on the momentum RHS",this,frame_num);
    write_scalar_attr(fh,"/fields/fp_total_z","fp_total is the actual particle-feedback acceleration used in the momentum RHS; z component cached during compute_ustar/acceleration, not reconstructed at output time","acceleration","positive z acceleration on the momentum RHS",this,frame_num);
    write_scalar_attr(fh,"/fields/fp_stressdiv_x","debug-only stress-divergence reconstruction div(sigma_total - sigma_fluid_pure), x component; not used as formal fp_total","acceleration","debug reconstruction",this,frame_num);
    write_scalar_attr(fh,"/fields/fp_stressdiv_y","debug-only stress-divergence reconstruction div(sigma_total - sigma_fluid_pure), y component; not used as formal fp_total","acceleration","debug reconstruction",this,frame_num);
    write_scalar_attr(fh,"/fields/fp_stressdiv_z","debug-only stress-divergence reconstruction div(sigma_total - sigma_fluid_pure), z component; not used as formal fp_total","acceleration","debug reconstruction",this,frame_num);
    write_scalar_attr(fh,"/fields/u_x","x component of carrier/mixture velocity stored on the solver cell-centered grid","velocity","positive x velocity",this,frame_num);
    write_scalar_attr(fh,"/fields/u_y","y component of carrier/mixture velocity stored on the solver cell-centered grid","velocity","positive y velocity",this,frame_num);
    write_scalar_attr(fh,"/fields/u_z","z component of carrier/mixture velocity stored on the solver cell-centered grid","velocity","positive z velocity",this,frame_num);
    write_scalar_attr(fh,"/fields/chi_f","1 for carrier-fluid cells and 0 for particle/solid cells; chi_f + chi_p = 1","dimensionless","mask",this,frame_num);
    write_scalar_attr(fh,"/fields/chi_p","1 for particle/solid cells and 0 for carrier-fluid cells; chi_f + chi_p = 1","dimensionless","mask",this,frame_num);
    write_scalar_attr(fh,"/fields/sigma_total_xx","total solver stress tensor component averaged from face stresses; mixed stress, not pure elastic sigma_e","stress-like solver units","diagnostic total stress",this,frame_num);
    write_attr(fh,"/diagnostics/K_full","budget_domain","full-domain one-fluid");
    write_attr(fh,"/diagnostics/K_fluid","budget_domain","carrier-fluid masked");
    write_attr(fh,"/diagnostics/dK_full_dt","definition","NaN online because future frame is unavailable; compute from K_full time series in post-processing");
    write_attr(fh,"/diagnostics/dK_fluid_dt","definition","NaN online because future frame is unavailable; compute from K_fluid time series in post-processing");
    write_attr(fh,"/diagnostics/residual_full","definition","dK_full_dt + eps_full - Wp_full; NaN online until dK_full_dt is available");
    write_attr(fh,"/diagnostics/residual_fluid","definition","dK_fluid_dt + eps_fluid - Wp_fluid; NaN online until dK_fluid_dt is available");
    write_attr(fh,"/diagnostics/ubar_full_x","average_definition","full-domain spatial average");
    write_attr(fh,"/diagnostics/ubar_fluid_x","average_definition","chi_f-weighted carrier-fluid spatial average");
    write_attr(fh,"/particles","particle_power_support","forcing/particle-interior nearest-center assigned using ppf3d_full and cached fp_rhs_total support");
    write_attr(fh,"/particles","deformation_available","0");
    write_attr(fh,"/particles","deformation_definition","unavailable in this solver output path; aspect ratio and principal axis are NaN rather than zero-filled");
    write_attr(fh,"/particles","alpha_near","1.5");
    write_attr(fh,"/particles","g0","1e-6");
    write_attr(fh,"/particles","p","1");
    write_attr(fh,"/particles","cluster_threshold","3");
    write_attr(fh,"/particles","periodic_distance_used","true");
    fclose(fh);
}
