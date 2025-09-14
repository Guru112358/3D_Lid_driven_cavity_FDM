#include "params.cpp"
#include "datastruct.cpp"
#include "preproc.cpp"
#include "postproc.cpp"
#include "discequations.cpp"



int main(int argc, char* argv[])
{

double residual=1;

int count=0;
bool loop_switch=true;

//init with junk values 
domain dom(0,0,0,0.0,0.0,0.0);

simparam sim(0.0,0,0.0,0.0,0.0,0.0,0.0);


read_params(argv[1],sim,dom);



flow_variables fv(dom.nx,dom.ny,dom.nz);

init_variables(fv.u,fv.v,fv.w,fv.p,fv.ustar,fv.vstar,fv.wstar,fv.pstar,dom);

//grid corner values
dmatrix uc(dom.nx,dom.ny,dom.nz);
dmatrix vc(dom.nx,dom.ny,dom.nz);
dmatrix wc(dom.nx,dom.ny,dom.nz);
dmatrix pc(dom.nx,dom.ny,dom.nz);


while(loop_switch)
{

 	solve_x_mom(fv.u,fv.ustar,fv.v,fv.w,fv.pstar,dom,sim);
	solve_y_mom(fv.v,fv.vstar,fv.u,fv.w,fv.pstar,dom,sim);
	solve_z_mom(fv.w,fv.wstar,fv.u,fv.v,fv.pstar,dom,sim);

	solve_pressure_correction_equation(fv.ustar,fv.vstar, fv.wstar,fv.p_prime,sim.pressure_iters,dom,sim);

	correct_p(fv.p,fv.pstar,fv.p_prime,dom,sim);

	correct_u(fv.unp1,fv.p_prime,fv.ustar,dom,sim);
	correct_v(fv.vnp1,fv.p_prime,fv.vstar,dom,sim);
	correct_w(fv.wnp1,fv.p_prime,fv.wstar,dom,sim);
	
	apply_bc(fv.unp1,fv.vnp1,fv.wnp1,fv.p,dom);
	

	swap_variables(fv.u,fv.v,fv.w, fv.p,fv.unp1,fv.vnp1,fv.wnp1,fv.pstar,dom);

	residual=compute_residual(fv.ustar,fv.vstar,fv.wstar,dom);
	
	if((residual<sim.tol&&(count!=0)))	
	{
	compute_collocated_values(fv.u,fv.v,fv.w,fv.p,uc,vc,wc,pc,dom);
	write_file_vtk_ascii(uc,vc,wc,pc,dom.nx,dom.ny,dom.nz,dom.dx,dom.dy,dom.dz);
	std::cout<<"converged to a tolerance of: "<<sim.tol<<" in "	<<count<<" Iterations"<<"\n";
	loop_switch=false;
	
	}
	
	
	if(count%sim.print_interval==0)
	{
	compute_collocated_values(fv.u,fv.v,fv.w,fv.p,uc,vc,wc,pc,dom);
	write_file_vtk_ascii(uc,vc,wc,pc,dom.nx,dom.ny,dom.nz,dom.dx,dom.dy,dom.dz);
	std::cout<<"|| iteration is: "<<count<<", continuity residual : " <<residual<<" ||"<<"\n";
	
	}
	
	if(count==sim.max_iters)
	{
	std::cout<<"Exiting, Maximum iterations reached: "<<sim.max_iters<<"\n";
	loop_switch=false;
	}
	
	
	count+=1;	

}

return 0;

}
