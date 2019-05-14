#include "./lbm/Domain.h"
#include <time.h>

struct myUserData
{
    std::ofstream oss_ss;
    double g;
    double nu;
    double R;
    double rhos;
};

void Report(LBM::Domain &dom, void *UD)
{
    myUserData &dat = (*static_cast<myUserData *> (UD));
    // size_t nx = dom.Ndim(0);
    // size_t ny = dom.Ndim(1);
    // size_t nz = dom.Ndim(2);
    // double dx = dom.dx;
    if(dom.Time <1e-6)
    {
        String fs;
        fs.Printf("%s.out","settling_hex_10.0_100.0");
        
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss<<Util::_10_6<<"Time"<<Util::_8s<<"Re"<<Util::_8s<<"Vx"<<Util::_8s<<"Vy"<<Util::_8s<<"X"<<Util::_8s<<"Y\n";
    }else{
        double v = dom.Particles.back().V(1);
        dat.oss_ss<<Util::_10_6<<dom.Time<<Util::_8s<<v*dat.R*2.0/dat.nu<<Util::_8s<<dom.Particles.back().V(0)<<Util::_8s<<dom.Particles.back().V(1)<<Util::_8s<<dom.Particles.back().X(0)<<Util::_8s<<dom.Particles.back().X(1)<<std::endl;
        
    }
}


double random(double a, double b)
{
    double rn = std::rand()/((double) RAND_MAX);
    return (b-a)*rn+a;
}
  
int main (int argc, char **argv) try
{
    
    std::srand((unsigned)time(NULL));    
    size_t Nproc = 12;
    double nu = 0.01;
    double ratio = 10.0;
    int pnx = 6;
    int pny = 70;
    double R = 5.0;
    double RR = ratio*R;
    bool initfromfile = false;
    char *h5name = NULL;
    if(argc>=2) Nproc = atoi(argv[1]);
    if(argc>=3) ratio = atof(argv[2]);
    if(argc>=4)
    {
        initfromfile = true;
        h5name = argv[3];
    }

    size_t nx = 600;
    size_t ny = 7000;
    size_t nz = 1;
    double dx = 1.0;
    double dt = 1.0;
    double Ga = 100.0;
    double rho = 1.0;
    double rhos = 1.1;
    double gy = Ga*Ga*nu*nu/((8*R*R*R)*(rhos/rho-1));
    std::cout<<"gy = "<<gy<<std::endl;

    std::cout<<"R = "<<R<<std::endl;
    //nu = 1.0/30.0;
    std::cout<<nx<<" "<<ny<<" "<<nz<<std::endl;
    LBM::Domain dom(D2Q9,MRT, nu, iVec3_t(nx,ny,nz),dx,dt);
    myUserData my_dat;
    dom.UserData = &my_dat;
    my_dat.nu = nu;
    my_dat.g = gy;
    my_dat.R = R;
    Vec3_t g0(0.0,0.0,0.0);
    dom.Nproc = Nproc;       

    //initial
    
    my_dat.rhos = rhos;
    
    
    Vec3_t pos(0.0,0.0,0.0);
    Vec3_t pos1(0.0,0.0,0.0);
    Vec3_t dxp(0.0,0.0,0.0);
    Vec3_t v(0.0,0.0,0.0);
    Vec3_t w(0.0,0.0,0.0);
    //DEM
    dom.dtdem = 0.01*dt;
    int pnum = 0;
    //fixed
    for(int j=0; j<pny; ++j)
    {
        int kx = 0;
        if(j%2!=0) kx = -1;
        for(int i = kx; i<pnx; ++i)
        {
            dxp = (2.0*i+j%2+1)*RR,(std::sqrt(3)*j+1)*RR,0.0;
            pos1 = pos+dxp;
            dom.Particles.push_back(DEM::Disk(pnum, pos1, v, w, rhos, 0.8*RR, dom.dtdem));
            dom.Particles[pnum].FixVeloc();
            
            pnum++;
        }
    }
    //move
    pos = nx/2-1,10.0,0;
    Vec3_t dxr(random(-0.1,0.1),random(-0.1,0.1));
    dom.Particles.push_back(DEM::Disk(-pnum, pos+dxr, v, w, rhos, R, dom.dtdem));
    pnum++;   

    std::cout<<"Particles number = "<<dom.Particles.size()<<std::endl;
    
    for(int ip=0; ip<(int) dom.Particles.size(); ++ip)
    {
        dom.Particles[ip].Ff = 0.0, M_PI*R*R*rhos*gy, 0.0;
        dom.Particles[ip].Kn = 20;
        dom.Particles[ip].Gn = 0.8;
        dom.Particles[ip].Kt = 2.5;
        // dom.Particles[ip].Mu = 0.0;
        // dom.Particles[ip].Eta = 0.0;
        // dom.Particles[ip].Beta = 0.0;
        dom.Particles[ip].Rh = dom.Particles[ip].R;
        // dom.Particles[ip].FixVeloc();

    }

    dom.Particles.back().Rh = 0.8*dom.Particles.back().R;
    
    for(size_t ix=0; ix<nx; ix++)
    {
        dom.IsSolid[ix][0][0] = true;
        dom.IsSolid[ix][ny-1][0] = true;
    }
    // for(size_t iy=0; iy<ny; iy++)
    // {
    //     dom.IsSolid[0][iy][0] = true;
    //     dom.IsSolid[nx-1][iy][0] = true;
    // }

    Vec3_t v0(0.0,0.0,0.0);
    dom.IsF = true;
    if(initfromfile)
    {
        dom.InitialFromH5(h5name,g0);

    }else{
        dom.Initial(rho,v0,g0);

    }


    double Tf = 1e6;
    
    double dtout = 1e3;
    dom.Box = 0.0,(double) nx-1, 0.0;
    dom.modexy = 0;
    //solving
    dom.Solve( Tf, dtout, "/home/Staff/uqyche38/macondo/settling_hex_10.0_100.0/test_fine_settling", NULL, Report);
    
    return 0;
}MECHSYS_CATCH
