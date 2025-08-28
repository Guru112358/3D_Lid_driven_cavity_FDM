

struct domain {
    int nx, ny, nz;       
    double lx, ly, lz;    
    double dx, dy, dz;    

    domain(int nx_, int ny_, int nz_, double lx_, double ly_, double lz_) {
        nx = nx_;
        ny = ny_;
        nz = nz_;
        lx = lx_;
        ly = ly_;
        lz = lz_;
        dx = lx / nx;
        dy = ly / ny;
        dz = lz / nz;
    }
};



struct simparam
{
    double dt;             
    double Re;            
    double tol;            
    int pressure_iters;    
    int print_interval;    
    long max_iters; 
    double urf_p;       

    simparam(double dt_,
             double Re_,
             double tol_,
             int pressure_iters_,
             int print_interval_,
             long max_iters_,
             double urf_p_)
        : dt(dt_),
          Re(Re_),
          tol(tol_),
          pressure_iters(pressure_iters_),
          print_interval(print_interval_),
          max_iters(max_iters_),
          urf_p(urf_p_)
    {}
};




struct flow_variables
{

dmatrix u,unp1,ustar;
dmatrix v,vnp1,vstar;
dmatrix w,wnp1,wstar;
dmatrix p,pnp1,p_prime,pstar;


flow_variables(int nx,int ny,int nz):u(nx,ny+1,nz+1),unp1(nx,ny+1,nz+1),ustar(nx,ny+1,nz+1),v(nx+1,ny,nz+1),vnp1(nx+1,ny,nz+1),vstar(nx+1,ny,nz+1),w(nx+1,ny+1,nz),wnp1(nx+1,ny+1,nz),wstar(nx+1,ny+1,nz),p(nx+1,ny+1,nz+1)
,pnp1(nx+1,ny+1,nz+1),p_prime(nx+1,ny+1,nz+1),pstar(nx+1,ny+1,nz+1) {} 

};
