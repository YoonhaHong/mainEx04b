#include <iostream>
#include <vector>
#include <string>
#define pi 3.141592
using namespace std;
class Pseudo:public TVector3
{
	public:
		Pseudo(float pt=0, float eta=0, float phi=0)
		{
			this->SetPtEtaPhi(pt,eta,phi);
		}
};

inline float getjt(Pseudo jet, Pseudo part) {return (jet.Cross(part)).Mag()/jet.Mag();}

inline float getz(Pseudo jet, Pseudo part) {return jet*part/jet.Mag2();}

float getr(float eta1, float eta2, float phi1, float phi2){ 
	float deta=eta1-eta2;
	float dphi=phi1-phi2;
	if(dphi<-pi) dphi+=2*pi;
	else if(dphi>pi) dphi-=2*pi;
	return sqrt(deta*deta+dphi*dphi);
}
void jetparticle2(int Rindex=5, TString inputpath="./allfull_MIF.root") //5 -> R=0.7
{
    float *p_pt=new float[5000];
	float *p_eta=new float[5000];
	float *p_phi=new float[5000];
	bool *p_chg=new bool[5000];
	int *p_id=new int[5000];
    int *p_status=new int[5000];

	float *ijet_pt=new float[20];
	float *ijet_eta=new float[20];
	float *ijet_phi=new float[20];

	int partnum;
	int ijetnum=0;
	int nijet=0,nqijet=0,ngijet=0;

    TH2F *qpt = new TH2F("quarkjet_particle_pt","quarkjet_particle_pt",1300,10,140,1300,10,140);
    TH2F *qeta = new TH2F("quarkjet_particle_eta","quarkjet_particle_eta",60,-3,3,80,-4,4);
    TH2F *qphi = new TH2F("quarkjet_particle_phi","quarkjet_particle_phi",60,-pi,pi,60,-pi,pi);
    TH1F *qdr = new TH1F("quarkdr","quarkdr",70,0.,0.7);

    TH2F *gpt = new TH2F("gluonjet_particle_pt","gluonjet_particle_pt",1300,10,140,1300,10,140);
    TH2F *geta = new TH2F("gluonjet_particle_eta","gluonjet_particle_eta",60,-3,3,80,-4,4);
    TH2F *gphi = new TH2F("gluonjet_particle_phi","gluonjet_particle_phi",60,-pi,pi,60,-pi,pi);
    TH1F *gdr = new TH1F("gluondr","gluondr",70,0.,0.7);

	TFile *f = new TFile(inputpath, "read");
	TTree *T = (TTree*)f ->Get("T");
	T->SetBranchAddress("p_pt",p_pt);
	T->SetBranchAddress("p_id",p_id);
    T->SetBranchAddress("p_status",p_status);
	T->SetBranchAddress("np",&partnum);
	T->SetBranchAddress("p_eta",p_eta);
	T->SetBranchAddress("p_phi",p_phi);
	T->SetBranchAddress("p_chg",p_chg);

    TString sijet;
	sijet.Form("ijetR%02d",Rindex+2);
	T->SetBranchAddress(sijet+"_pt",ijet_pt);
    T->SetBranchAddress(sijet+"_eta",ijet_eta);
    T->SetBranchAddress(sijet+"_phi",ijet_phi);
    T->SetBranchAddress("n"+sijet,&ijetnum);

	
    
	float jet_R_size = 0.1*(2.0+(float)Rindex); //R=0.7(Rindex=5) eta=-0.2~0.2
	float leading_pt; 
    int leading_index;
    int tagged_id;
    float dr;

	const int entrysum=T->GetEntries();

	for (int ien=0; ien<entrysum; ien++){
		T->GetEntry(ien);

		if ( ijetnum<1 ) continue;
        
		for (int ijet=0; ijet<ijetnum; ijet++){
            if(ijet_pt[ijet]>140 or ijet_pt[ijet]<10) continue;
            if(fabs(ijet_eta[ijet])>3) continue;
            leading_pt=0;
			for (int part=0; part<partnum; part++){
                
				dr = getr(ijet_eta[ijet],p_eta[part],ijet_phi[ijet],p_phi[part]);
                if ( dr > 0.7) continue; //for _07finding
				if( p_pt[part]>leading_pt) {
					leading_pt = p_pt[part];
					leading_index = part; 
				}
			}//part
            nijet++;

            tagged_id = fabs( p_id[leading_index] );
            if( tagged_id==1 or tagged_id==2 or tagged_id==3){
                nqijet++;
               
                qpt->Fill(ijet_pt[ijet], p_pt[leading_index]);
                qeta->Fill(ijet_eta[ijet], p_eta[leading_index]);
                qphi->Fill(ijet_phi[ijet], p_phi[leading_index]);
                qdr->Fill( getr(ijet_eta[ijet],p_eta[leading_index],ijet_phi[ijet],p_phi[leading_index]) );
            }
            else{
                ngijet++;

                gpt->Fill(ijet_pt[ijet], p_pt[leading_index]);
                geta->Fill(ijet_eta[ijet], p_eta[leading_index]);
                gphi->Fill(ijet_phi[ijet], p_phi[leading_index]);
                gdr->Fill( getr(ijet_eta[ijet],p_eta[leading_index],ijet_phi[ijet],p_phi[leading_index]) );
            }
            
		}//ijet
    }//ien
    cout<<"quark tagged: "<<nqijet<<endl
        <<"gluon tagged(not quark): "<<ngijet<<endl
        <<"all jets in condition: "<<nijet<<endl;
    
    //return;

	TString outpath="./jetparticle2_all_MIF_";
   /*
    if(inputpath.Contains("MPI")) outpath="./MPI/jetparticle";
    else if(inputpath.Contains("ISR")) outpath="./ISR/jetparticle";
    else if(inputpath.Contains("FSR")) outpath="./FSR/jetparticle";

    if(mode==1) outpath.Append("gg2qqbar");
    else outpath.Append("gg2gg");

    if(inputpath.Contains("MIF")) outpath.Append("MIF");
    else if(inputpath.Contains("MPI")) outpath.Append("MPI");
    else if(inputpath.Contains("ISR")) outpath.Append("ISR");
    else if(inputpath.Contains("FSR")) outpath.Append("FSR");
*/
    TString Rsize;
	Rsize.Form("R%02d",Rindex+2);
    outpath.Append(Rsize);
    outpath.Append(".root");
    TFile *outfile = new TFile(outpath, "RECREATE");
    
    qpt->Write();
    qeta->Write();
    qphi->Write();
    qdr->Write();

    gpt->Write();
    geta->Write();
    gphi->Write();
    gdr->Write();

    outfile->Write();
    outfile->Close();
    return;
}
