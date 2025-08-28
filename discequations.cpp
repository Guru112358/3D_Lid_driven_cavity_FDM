void set_zero(dmatrix &A);
void init_variables(dmatrix &u,dmatrix&v, dmatrix &w ,dmatrix &p,dmatrix &ustar,dmatrix&vstar, dmatrix &wstar, dmatrix &pstar,domain&dom);
void solve_x_mom(dmatrix&u,dmatrix &unp1,dmatrix &v,dmatrix &w,dmatrix &p,domain&dom,simparam&sim);
void solve_y_mom(dmatrix&v,dmatrix &vnp1,dmatrix &u,dmatrix &w,dmatrix &p,domain&dom,simparam&sim);
void solve_z_mom(dmatrix&w,dmatrix &wnp1,dmatrix &u,dmatrix &v,dmatrix &p,domain&dom,simparam&sim);
void apply_bc(dmatrix&u,dmatrix &v, dmatrix &w, dmatrix &p,domain&dom);
double compute_residual(dmatrix&u,dmatrix &v,dmatrix&w,domain&dom);
void compute_collocated_values(dmatrix &u,dmatrix &v , dmatrix &w ,dmatrix &p,dmatrix &uc,dmatrix &vc, dmatrix &wc ,dmatrix &pc,domain&dom);
void swap_variables(dmatrix &u,dmatrix &v ,dmatrix &p,dmatrix &unp1,dmatrix &vnp1,dmatrix &pnp1,dmatrix &ustar,dmatrix &vstar,dmatrix &pstar,domain&dom);
void solve_pressure_correction_equation(dmatrix&u,dmatrix &v, dmatrix &w,dmatrix&p_prime,int niter,domain&dom,simparam&sim);
void correct_u(dmatrix &unp1,dmatrix &pprime,dmatrix &ustar,domain&dom,simparam&sim);
void correct_v(dmatrix &vnp1,dmatrix &pprime,dmatrix &vstar,domain&dom,simparam&sim);
void correct_w(dmatrix &wnp1,dmatrix &pprime,dmatrix &wstar,domain&dom,simparam&sim);
void correct_p(dmatrix &p,dmatrix &pstar,dmatrix&p_prime,domain&dom,simparam&sim);


void init_variables(dmatrix &u,dmatrix&v, dmatrix &w ,dmatrix &p,dmatrix &ustar,dmatrix&vstar, dmatrix &wstar, dmatrix &pstar,domain&dom)
{


for(int i=0;i<=(dom.nx-1);i++)
{
	for(int j=0;j<=dom.ny;j++)
	{	
	   for(int k=0;k<=dom.nz;k++)
    	{
    		ustar(i,j,k)=0;
    		u(i,j,k)=0.0;
		
    	}
	}
}


	
for(int i=0;i<=dom.nx;i++)
{
	for(int j=0;j<=(dom.ny-1);j++)
	{	
	  for(int k=0;k<=dom.nz;k++)
	    {
		    v(i,j,k)=0.0;
	    	    vstar(i,j,k)=0;
		
    	}
	}
}


for(int i=0;i<=dom.nx;i++)
{
	for(int j=0;j<=(dom.ny);j++)
	{	
	  for(int k=0;k<=(dom.nz-1);k++)
	    {
		    w(i,j,k)=0.0;
	    	    wstar(i,j,k)=0;
		
    	}
	}
}


for(int i=0;i<=dom.nx;i++)
{
	for(int j=0;j<=dom.ny;j++)
	{	 
	    for(int k=0;k<=dom.nz;k++)
	    {
		    p(i,j,k)=0;
		    pstar(i,j,k)=1;
		}
	}
}

}


void solve_x_mom(dmatrix&u,dmatrix &unp1,dmatrix &v,dmatrix &w,dmatrix &p,domain&dom,simparam &sim)
{

double Ax,Ay,Az,Dx,Dy,Dz;
double gradp_x;


#pragma omp parallel for collapse(3) schedule(static) private(gradp_x,Ax,Ay,Az,Dx,Dy,Dz)

for(int i=1;i<=(dom.nx-2);i++)
{
	for(int j=1;j<=(dom.ny-1);j++)
	{
		for(int k=1;k<=(dom.nz-1);k++)
		{

	    gradp_x=p(i+1,j,k)-p(i,j,k);

		Ax=(u(i+1,j,k)*u(i+1,j,k) - u(i-1,j,k)*u(i-1,j,k)) / (2.0*dom.dx);
		
		Ay=0.25 *(((u(i,j,k) + u(i,j+1,k)) * (v(i,j,k) + v(i+1,j,k)) - (u(i,j,k) + u(i,j-1,k)) * (v(i+1,j-1,k) + v(i,j-1,k)) ) / dom.dy);

		Az= 0.25*(( ((u(i,j,k)+u(i,j,k+1))*(w(i,j,k)+w(i+1,j,k)) ) - ((u(i,j,k)+u(i,j,k-1))*(w(i,j,k-1)+w(i+1,j,k-1))) )/dom.dz);

		Dx=(u(i+1,j,k)+u(i-1,j,k)-2*u(i,j,k))/(dom.dx*dom.dx);
		
		Dy=(u(i,j+1,k)+u(i,j-1,k)-2*u(i,j,k))/(dom.dy*dom.dy);

		Dz=(u(i,j,k+1)+u(i,j,k-1)-2*u(i,j,k))/(dom.dz*dom.dz);
	
		unp1(i,j,k)=u(i,j,k) - (Ax+Ay+Az)*sim.dt + (sim.dt/sim.Re)*(Dx+Dy+Dz) - (sim.dt/dom.dx)*gradp_x;


		}	
	
	}

}

}


void solve_y_mom(dmatrix&v,dmatrix &vnp1,dmatrix &u,dmatrix &w,dmatrix &p,domain&dom,simparam&sim)
{

double Ax,Ay,Az,Dx,Dy,Dz;
double gradp_y;

#pragma omp parallel for collapse(3) schedule(static) private(gradp_y,Ax,Ay,Az,Dx,Dy,Dz)

for (int i = 1; i <= dom.nx - 1; i++)
{
   for (int j = 1; j <= dom.ny - 2; j++)
	 {
		for(int k=1; k<=(dom.nz-1); k++)
		{

		gradp_y=p(i,j+1,k)-p(i,j,k);

		Ax=0.25*(((u(i,j,k) + u(i,j+1,k)) * (v(i,j,k) + v(i+1,j,k)) - (u(i-1,j,k) + u(i-1,j+1,k)) * (v(i,j,k) + v(i-1,j,k)))/ dom.dx);

		Ay=(v(i,j+1,k) * v(i,j+1,k) - v(i,j-1,k) * v(i,j-1,k)) / (2.0 * dom.dy);

		Az= 0.25*(( ((v(i,j,k)+v(i,j,k+1))*(w(i,j,k)+w(i,j+1,k)) )    -((v(i,j,k)+v(i,j,k-1))*(w(i,j,k-1)+w(i,j+1,k-1))) )/dom.dz);

		Dx=(v(i+1,j,k)+v(i-1,j,k)-2*v(i,j,k))/(dom.dx*dom.dx);
		
		Dy=(v(i,j+1,k)+v(i,j-1,k)-2*v(i,j,k))/(dom.dy*dom.dy);	

		Dz=(v(i,j,k+1)+v(i,j,k-1)-2*v(i,j,k))/(dom.dz*dom.dz);
		
		vnp1(i,j,k)=v(i,j,k) - (Ax+Ay+Az)*sim.dt + (sim.dt/sim.Re)*(Dx+Dy+Dz) - (sim.dt/dom.dy)*gradp_y;

		}
   }
}

}


void solve_z_mom(dmatrix&w,dmatrix &wnp1,dmatrix &u,dmatrix &v,dmatrix &p,domain&dom,simparam&sim)
{

double Ax,Ay,Az,Dx,Dy,Dz;
double gradp_z;

#pragma omp parallel for collapse(3) schedule(static) private(gradp_z,Ax,Ay,Az,Dx,Dy,Dz)

for (int i = 1; i <= dom.nx - 1; i++)
{
   for (int j = 1; j <= dom.ny - 1; j++)
	 {
		for (int k = 1; k <= dom.nz - 2; k++)
		{

		gradp_z=p(i,j,k+1)-p(i,j,k);

		Ax=0.25*(((u(i,j,k) + u(i,j,k+1)) * (w(i,j,k) + w(i+1,j,k)) - (u(i-1,j,k) + u(i-1,j,k+1)) * (w(i,j,k) + w(i-1,j,k)))/ dom.dx);

		Ay=0.25 *(((w(i,j,k) + w(i,j+1,k)) * (v(i,j,k) + v(i,j,k+1)) - (w(i,j,k) + w(i,j-1,k)) * (v(i,j-1,k) + v(i,j-1,k+1)) ) / dom.dy );

		Az=(w(i,j,k+1) * w(i,j,k+1) - w(i,j,k-1) * w(i,j,k-1)) / (2.0 * dom.dz);

		Dx=(w(i+1,j,k)+w(i-1,j,k)-2*w(i,j,k))/(dom.dx*dom.dx);
		
		Dy=(w(i,j+1,k)+w(i,j-1,k)-2*w(i,j,k))/(dom.dy*dom.dy);	

		Dz=(w(i,j,k+1)+w(i,j,k-1)-2*w(i,j,k))/(dom.dz*dom.dz);

		wnp1(i,j,k)=w(i,j,k)-(Ax+Ay+Az)*sim.dt +(sim.dt/sim.Re)*(Dx+Dy+Dz) -(sim.dt/dom.dz)*gradp_z;

		
   }
}

}

}


void solve_pressure_correction_equation(dmatrix&u,dmatrix &v, dmatrix &w,dmatrix&p_prime,int niter,domain&dom,simparam&sim)
{


double a,b,c,d,e,rhs;

a=2.0*( ((sim.dt)/(dom.dx*dom.dx))  + ((sim.dt)/(dom.dy*dom.dy))  +  ((sim.dt)/(dom.dz*dom.dz)));
b=sim.dt/(dom.dx*dom.dx);
c=sim.dt/(dom.dy*dom.dy);
e=sim.dt/(dom.dz*dom.dz);

double norm;

set_zero(p_prime);

for(int iter_count=0;iter_count<niter;iter_count++)
{

	d=0;
	rhs=0;
	
	
	#pragma omp for collapse(3) schedule(static) private(d,rhs) 
	

	for(int i = 1; i <= (dom.nx - 1); i++)
	{
	  for(int j = 1; j <= (dom.ny - 1); j++)
		{
		for(int k=1;k <= (dom.nz-1); k++)
			{

			d= (u(i,j,k) - u(i-1,j,k)) / dom.dx + (v(i,j,k) - v(i,j-1,k)) / dom.dy  + (w(i,j,k)-w(i,j,k-1))/dom.dz;  //mass balance source term;

			rhs=(b*(p_prime(i+1,j,k)+p_prime(i-1,j,k))+c*(p_prime(i,j+1,k)+p_prime(i,j-1,k)) + e*(p_prime(i,j,k+1)+p_prime(i,j,k-1))-d)/a;

			p_prime(i,j,k)=rhs;


			
		}
	}
}

}
}

void correct_u(dmatrix &unp1,dmatrix &pprime,dmatrix &ustar,domain&dom,simparam&sim)
{

double gradp_x;



for(int i=1;i<=(dom.nx-2);i++)
{
	for(int j=1;j<=(dom.ny-1);j++)
	{	
		for(int k=1;k<=(dom.nz-1);k++)
		{
		
    		gradp_x=pprime(i+1,j,k)-pprime(i,j,k);
	  		unp1(i,j,k)=ustar(i,j,k)-(sim.dt/dom.dx)*gradp_x;
		}

	}
}

}

void correct_v(dmatrix &vnp1,dmatrix &pprime,dmatrix &vstar,domain&dom,simparam&sim)
{

double gradp_y;



for (int i = 1; i <= dom.nx - 1; i++)
	{
   	for (int j = 1; j <= dom.ny - 2; j++)
	 {
		for(int k=1; k <= (dom.nz-1); k++)
			{
	 
     		 	 gradp_y=pprime(i,j+1,k)-pprime(i,j,k);
	 		  vnp1(i,j,k)=vstar(i,j,k)-(sim.dt/dom.dy)*gradp_y;

			}
		}
	}
}



void correct_w(dmatrix &wnp1,dmatrix &pprime,dmatrix &wstar,domain&dom,simparam&sim)
{

double gradp_z;



for (int i = 1; i <= dom.nx - 1; i++)
	{
   	for (int j = 1; j <= dom.ny - 1; j++)
	 {
		for(int k=1; k <= (dom.nz-2); k++)
			{
	 
     			  gradp_z=pprime(i,j,k+1)-pprime(i,j,k);
	 		  wnp1(i,j,k)=wstar(i,j,k)-(sim.dt/dom.dz)*gradp_z;

			}
		}
	}
}



void correct_p(dmatrix &p,dmatrix &pstar,dmatrix&p_prime,domain&dom,simparam&sim)
{


for(int i = 1; i <= (dom.nx - 1); i++)
{
	for(int j = 1; j <= (dom.ny - 1); j++)
	{
		for(int k=1;k <= (dom.nz-1); k++)
		{
			p(i,j,k)=pstar(i,j,k)+sim.urf_p*p_prime(i,j,k);
		}
	}
}

}






void apply_bc(dmatrix&u,dmatrix &v, dmatrix &w, dmatrix &p,domain&dom)
{

//################# front plane ##################

//u velocity 


for(int i=0;i<=dom.nx-1;i++)
	for(int j=0;j<=dom.ny;j++)
{
	u(i,j,dom.nz)=-u(i,j,dom.nz-1);
}
//v velocity

for(int i=0;i<=dom.nx;i++)
	for(int j=0;j<=dom.ny-1;j++)
{
	v(i,j,dom.nz)=-v(i,j,dom.nz-1);
}

//w velocity

for(int i=0;i<=dom.nx;i++)
	for(int j=0;j<=dom.ny;j++)
{
	w(i,j,dom.nz-1)=0;
}

//pressure

for(int i=0;i<=dom.nx;i++)
	for(int j=0;j<=dom.ny;j++)
{
	p(i,j,dom.nz)=p(i,j,dom.nz-1);
}


//############ back plane ################

//u velocity 

for(int i=0;i<=dom.nx-1;i++)
	for(int j=0;j<=dom.ny;j++)
{
	u(i,j,0)=-u(i,j,1);
}
//v velocity

for(int i=0;i<=dom.nx;i++)
	for(int j=0;j<=dom.ny-1;j++)
{
	v(i,j,0)=-v(i,j,1);
}

//w velocity

for(int i=0;i<=dom.nx;i++)
	for(int j=0;j<=dom.ny;j++)
{
	w(i,j,0)=0;
}

//pressure

for(int i=0;i<=dom.nx;i++)
	for(int j=0;j<=dom.ny;j++)
{
	p(i,j,0)=p(i,j,1);
}
 



//################### bottom plane #####################


//u velocity 

for(int i=0;i<=dom.nx-1;i++)
	for(int k=0;k<=dom.nz;k++)
{
	u(i,0,k)=-u(i,1,k);
}
//v velocity

for(int i=0;i<=dom.nx;i++)
	for(int k=0;k<=dom.nz;k++)
{	

    v(i,0,k)=0;
}

//w velocity

for(int i=0;i<=dom.nx;i++)
	for(int k=0;k<=dom.nz-1;k++)
{
		w(i,0,k)=-w(i,1,k);

}

//pressure

for(int i=0;i<=dom.nx;i++)
	for(int k=0;k<=dom.nz;k++)
{
		p(i,0,k)=p(i,1,k);
}



// ############################# left plane ############################# 

//u velocity 

for(int j=0;j<=dom.ny;j++)
	for(int k=0;k<=dom.nz;k++)
{
	u(0,j,k)=0;
}
//v velocity

for(int j=0;j<=dom.ny-1;j++)
	for(int k=0;k<=dom.nz;k++)
{
	v(0,j,k)=-v(1,j,k);
}

//w velocity

for(int j=0;j<=dom.ny;j++)
	for(int k=0;k<=dom.nz-1;k++)
{
	w(0,j,k)=-w(1,j,k);
}

//pressure

for(int j=0;j<=dom.ny;j++)
	for(int k=0;k<=dom.nz;k++)
{
	p(0,j,k)=p(1,j,k);
}


// ############## right plane #######################

//u velocity 

for(int j=0;j<=dom.ny;j++)
	for(int k=0;k<=dom.nz;k++)
{
	u(dom.nx-1,j,k)=0;
}
//v velocity

for(int j=0;j<=dom.ny-1;j++)
	for(int k=0;k<=dom.nz;k++)
{
	v(dom.nx,j,k)=-v(dom.nx-1,j,k);
}

//w velocity

for(int j=0;j<=dom.ny;j++)
	for(int k=0;k<=dom.nz-1;k++)
{
	w(dom.nx,j,k)=-w(dom.nx-1,j,k);
}

//pressure

for(int j=0;j<=dom.ny;j++)
	for(int k=0;k<=dom.nz;k++)
{
	p(dom.nx,j,k)=p(dom.nx-1,j,k);
}






//////CRITICAL SECTION //////////////

//############## top plane #######################

//u velocity 

for(int i=0;i<=dom.nx-1;i++)
	for(int k=0;k<=dom.nz;k++)
{
	u(i,dom.ny,k)=2-u(i,dom.ny-1,k);
}
//v velocity

for(int i=0;i<=dom.nx;i++)
	for(int k=0;k<=dom.ny;k++)
{	

	
	v(i,dom.ny-1,k)=0;

}

//w velocity

for(int i=0;i<=dom.nx;i++)
	for(int k=0;k<=dom.nz-1;k++)
{
    w(i,dom.ny,k)=-w(i,dom.ny-1,k);


}

//pressure

for(int i=0;i<=dom.nx;i++)
	for(int k=0;k<=dom.nz;k++)
{
		p(i,dom.ny,k)=p(i,dom.ny-1,k);
}

}






double compute_residual(dmatrix&u,dmatrix &v,dmatrix&w,domain&dom)
{
double res=0;

double mb;

#pragma omp parallel for schedule(static) private(mb) reduction(+:res)

for (int i = 1; i <= (dom.nx - 1); i++)
{
	for (int j = 1; j <= (dom.ny - 1); j++)
	{
		for (int k = 1; k <= (dom.nz - 1); k++)
		{
			mb = (u(i,j,k) - u(i-1,j,k)) / dom.dx + (v(i,j,k) - v(i,j-1,k)) / dom.dy   +(w(i,j,k)-w(i,j,k-1))/dom.dz;
			res+= fabs(mb);
		
		}
	}
}

return res;

}



void compute_collocated_values(dmatrix &u,dmatrix &v , dmatrix &w ,dmatrix &p,dmatrix &uc,dmatrix &vc, dmatrix &wc ,dmatrix &pc,domain&dom)
{

#pragma omp for collapse(3) 

for (int i = 0; i <= (dom.nx - 1); i++)
{
	for (int j = 0; j <= (dom.ny - 1); j++)
	{
		
		for (int k = 0; k <= (dom.nz - 1); k++)
		{	
		
			uc(i,j,k) = 0.25 * (u(i,j,k) + u(i,j+1,k) + u(i,j,k+1) + u(i,j+1,k+1)  );

			vc(i,j,k) = 0.25 * (v(i,j,k) + v(i+1,j,k) + v(i,j,k+1) + v(i+1,j,k+1) );

          		wc(i,j,k)=  0.25* (w(i+1,j,k)+w(i,j+1,k)+w(i+1,j+1,k)+w(i,j,k));
				
			pc(i,j,k) =   (1.0/8) * (p(i,j,k) + p(i+1,j,k) + p(i,j+1,k) + p(i+1,j+1,k) +

			p(i,j,k+1) + p(i+1,j,k+1) + p(i,j+1,k+1) + p(i+1,j+1,k+1)) ;

		}
	}
}

}




void swap_variables(dmatrix &u,dmatrix &v , dmatrix &w, dmatrix &p,dmatrix &unp1,dmatrix &vnp1, dmatrix &wnp1,dmatrix &pstar,domain&dom)
{



for(int i=0;i<=(dom.nx-1);i++)
{
	for(int j=0;j<=dom.ny;j++)
	{	
	   for(int k=0;k<=dom.nz;k++)
    	{

			u(i,j,k)=unp1(i,j,k);
		}
	}
}


	
for(int i=0;i<=dom.nx;i++)
{
	for(int j=0;j<=(dom.ny-1);j++)
	{	
	  for(int k=0;k<=dom.nz;k++)
	    {
			v(i,j,k)=vnp1(i,j,k);
		}
	}
}


for(int i=0;i<=dom.nx;i++)
{
	for(int j=0;j<=(dom.ny);j++)
	{	
	  for(int k=0;k<=(dom.nz-1);k++)
	    {
			w(i,j,k)=wnp1(i,j,k);
		}
	}
}




for (int i = 0; i <= dom.nx; i++)
{
	for (int j = 0; j <= dom.ny; j++)
	{
		for(int k=0;k<=dom.nz;k++)
		{
			pstar(i,j,k)=p(i,j,k);
		}
	}
}

}



void set_zero(dmatrix &A)
{
#pragma omp parallel for collapse(3)
	for (int i = 0; i < A.dimension(0); i++)
{
	for (int j = 0; j <A.dimension(1); j++)
	{

		for(int k=0;k<A.dimension(2);k++)
		{
			A(i,j,k)=0;
		}
	}
}

}
