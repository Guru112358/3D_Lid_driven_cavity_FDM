void read_params(const std::string &filename,simparam&sim,domain&dom);





struct simparam
{
    double dt;             
    double Re;            
    double tol;            
    int pressure_iters;    
    int print_interval;    
    long max_iters; 
    double urf_p;       
    int nthreads;
    simparam(double dt_,
             double Re_,
             double tol_,
             int pressure_iters_,
             int print_interval_,
             long max_iters_,
             double urf_p_,int n_threads)
        : dt(dt_),
          Re(Re_),
          tol(tol_),
          pressure_iters(pressure_iters_),
          print_interval(print_interval_),
          max_iters(max_iters_),
          urf_p(urf_p_),
          nthreads(n_threads)
    {}
};