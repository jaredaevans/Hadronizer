// Generate LHC collisions and extract hadronization only

#define DO_ISOLATION        true
#define WRITE_MATH          false
#define WRITE_VERT          false

#define PI 3.14159265358979323846


#include "Pythia8/Pythia.h"
using namespace Pythia8;


/////////// START: Routines for Delta R /////////////
float sqr(float x)
{ return x*x; }
float getEta(float px, float py, float pz)
{
    float p = sqrt(sqr(px) + sqr(py) + sqr(pz));
    if (p != pz)
        return log((p + pz) / (p - pz)) / 2;
    else
        return 999.;
}
float getPhi(float px, float py)  // returns number between 0 and 2pi
{
    float phi;
    if (px == 0)
        phi = (py > 0 ? PI/2 : 3*PI/2);
    else {
        phi = atan2(py, px);
        if (phi < 0) phi += 2*PI;
    }
    return phi;
}
float phi02pi(float phi) // converts phi to be between 0 and 2pi
{ return (phi < 0 ? phi + 2*PI : phi); }

float dphi(float p1x, float p1y, float p2x, float p2y)
{
    float dphi = fabs(getPhi(p1x, p1y) - getPhi(p2x, p2y));
    if (dphi > PI)  dphi = 2*PI - dphi;
    return dphi;
}
float dR(float p1x, float p1y, float p1z, float p2x, float p2y, float p2z)
{ return sqrt(sqr(getEta(p1x, p1y, p1z) - getEta(p2x, p2y, p2z)) + sqr(dphi(p1x, p1y, p2x, p2y))); }

float dR(Particle P1, Particle P2)
{ return dR(P1.px(), P1.py(), P1.pz(), P2.px(), P2.py(), P2.pz()); }
/////////// END: Routines for Delta R /////////////

/////////// START: Routines for B and D ID /////////////
/*bool IDisB(const Particle P)
 {
 if ((P.idAbs()>500 && P.idAbs()<600) || (P.idAbs()>5000 && P.idAbs()<6000))
 return true;
 return false;
 }
 
 bool IDisD(const Particle P)
 {
 if ((P.idAbs()>400 && P.idAbs()<500) || (P.idAbs()>4000 && P.idAbs()<5000))
 return true;
 return false;
 }*/
bool IDisB(const Particle P)
{
    if (P.idAbs()==511 || P.idAbs()==521 || P.idAbs()==531|| (P.idAbs()>5000 && P.idAbs()<6000))
        return true;
    return false;
}

bool IDisD(const Particle P)
{
    if (P.idAbs()==411 || P.idAbs()==421  || P.idAbs()==431 || (P.idAbs()>4000 && P.idAbs()<5000))
        return true;
    return false;
}
/////////// END: Routines for B and D ID ////////////

/////////// START: Routines for I/O /////////////
void writeParticle(const Particle P,ofstream& o)
{
    o.precision(9);
#if WRITE_MATH
    o << "{" << P.id()
    <<  "," << P.e()
    <<  "," << P.px()
    <<  "," << P.py()
    <<  "," << P.pz()
    <<  "," << P.m()
    <<  "}";
#else
    o << P.id()
    <<  " " << P.e()
    <<  " " << P.px()
    <<  " " << P.py()
    <<  " " << P.pz()
#if WRITE_VERT
    <<  " " << P.m()
    <<  " " << P.xProd()
    <<  " " << P.yProd()
    <<  " " << P.zProd()
    <<  " " << P.tProd();
#else
    <<  " " << P.m();
#endif
#endif
    return;
}
void writeParticle(const Particle P,ostream& o)
{
    o.precision(9);
    o << "{" << P.id()
    <<  "," << P.e()
    <<  "," << P.px()
    <<  "," << P.py()
    <<  "," << P.pz()
    <<  "," << P.m()
    <<  "}";
    return;
}

void writeParticleList(const vector<Particle> Plist,ofstream& o)
{
    if(Plist.size()==0)
    {
        cout << "Zero size" << endl;
        return;
    }
#if WRITE_MATH
    o << "{";
    for(auto i=Plist.begin(); i<Plist.end()-1; i++)
    {
        writeParticle(*i,o);
        o << ",";
    }
    writeParticle(*(Plist.end()-1),o);
    o <<  "}";
#else
    for(auto i=Plist.begin(); i<Plist.end()-1; i++)
    {
        writeParticle(*i,o);
        o << endl;
    }
    writeParticle(*(Plist.end()-1),o);
    o <<  endl;
#endif
    return;
}
void writeParticleList(const vector<Particle> Plist,ostream& o)
{
    if(Plist.size()==0)
    {
        cout << "Zero size" << endl;
        return;
    }
    o << "{";
    for(auto i=Plist.begin(); i<Plist.end()-1; i++)
    {
        writeParticle(*i,o);
        o << ",";
    }
    writeParticle(*(Plist.end()-1),o);
    o <<  "}";
    return;
}

void WriteHadronizedParticles(const vector<Particle> Plist,const vector<Particle> Qlist,ofstream& o)
{
#if WRITE_MATH
    o << "{";
    writeParticleList(Plist,o);
    o << ",";
    writeParticleList(Qlist,o);
    o <<  "}";
#else
    o << "&";
    writeParticleList(Plist,o);
    o << ";";
    writeParticleList(Qlist,o);
#endif
    return;
}
void WriteHadronizedParticles(const vector<Particle> Plist,const vector<Particle> Qlist,ostream& o)
{
    o << "{";
    writeParticleList(Plist,o);
    o << ",";
    writeParticleList(Qlist,o);
    o <<  "}";
    return;
}
/////////// END: Routines for I/O /////////////

bool contains(vector<int> mylist, int element) {
    for (int i = 0; i<mylist.size(); ++i)
    {
        if (mylist[i] == element)
            return true;
    }
    return false;
}

//==========================================================================

int main(int argc,char *argv[]) {
    if(argc!=3)
    {
        cout << "Form:  ./{EXE} {mode} {rnd seed}" << endl << endl;
        return 1;
    }
    
    char command[100];
    int mode = atoi(argv[1]);
    int random_seed = atoi(argv[2]);
    
    int Nevents = 1 * mode;
    //int Nevents = mode;
    
    // Generator.
    Pythia pythia;
    //Event& event      = pythia.event;
    //ParticleData& pdt = pythia.particleData;
    
    string exptail = "";
    
    // LHC initialization
    pythia.readString("Beams:frameType = 1");
    pythia.readString("Beams:eCM =  13000"); // fixed to 13 TeV
    pythia.readString("HardQCD:all = on");
    
    string filetail = "";
    
    // Switch off all Z0 decays and then switch back on those to quarks.
    //pythia.readString("23:onMode = off");

    pythia.rndm.init(random_seed);
    
    pythia.init();

    char filename[100];
    sprintf(filename,"Hadronization.dat");
    ofstream outfileB(filename);
    
    int nNu = 0;
    int iAbort = 0;
    int nAbort=100000;
    
    vector<Particle> particleListx;
    vector<Particle> particleListy;
    
    int startval = 0;
    int endval = 0;
    bool writeAll = true;
    bool writeCur = true;
    
    
    vector<int> listofwritten;
    
    // Begin event loop. Generate event. Skip if error. List first one.
    for (int iEvent = 0; iEvent < Nevents; ++iEvent) {
        if (!pythia.next()) {
            if (++iAbort < nAbort)
            {
                iEvent--;
                continue;
            }
            cout << " Event generation aborted prematurely, owing to error!\n";
            break;
        }
        listofwritten.clear();
        
        for (int i = 0; i < pythia.event.size(); ++i)
        {
           Particle Q = pythia.event[i];
           if(Q.status() == -71)
           {
               if(!contains(listofwritten,i))
               {
                   for(int j = Q.daughter1(); j<= Q.daughter2(); ++j)
                       particleListy.push_back(pythia.event[j]);
                   Particle P = pythia.event[Q.daughter1()];
                   for(int j = P.mother1(); j<= P.mother2(); ++j)
                   {
                       particleListx.push_back(pythia.event[j]);
                       listofwritten.push_back(j);
                   }
                   WriteHadronizedParticles(particleListx,particleListy,outfileB);
                   particleListx.clear();
                   particleListy.clear();
               }
           }
        }
    }
    
    pythia.stat();
    
    double xsec = pythia.info.sigmaGen();
    double xsecerror = pythia.info.sigmaErr();
    
    outfileB.close();
    
    
    return 0;
}
