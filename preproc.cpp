void read_params(const std::string &filename,simparam&sim,domain&dom);




void read_params(const std::string &filename,simparam&sim,domain&dom) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open " << filename << "\n";
        std::exit(1);
}

file>>dom.lx;
file>>dom.ly;
file>>dom.lz;
file>>dom.nx;
file>>dom.ny;
file>>dom.nz;

dom.dx=dom.lx/dom.nx;
dom.dy=dom.ly/dom.ny;
dom.dz=dom.lz/dom.nz;

file>>sim.dt;
file>>sim.Re;
file>>sim.tol;
file>>sim.pressure_iters;
file>>sim.print_interval;
file>>sim.max_iters;
file>>sim.urf_p;
file.close();

}   
