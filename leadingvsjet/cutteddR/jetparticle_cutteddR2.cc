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
void jetparticle_cutteddR2(int Rindex=0, TString inputpath="../gg2ggfull_ISR.root") //5 -> R=0.7
{
    int mode=inputpath.Contains("qq"); //if quark jet mode =1 
    
    float ptcut[2] = {15, 20}; //jet pT>15, leading particle pT>20
    float *p_pt = new float[2000];
	float *p_eta = new float[2000];
	float *p_phi = new float[2000];
	bool *p_chg = new bool[2000];
	int *p_id = new int[2000];
    int *p_status = new int[2000];

	float *ijet_pt = new float[20];
	float *ijet_eta = new float[20];
	float *ijet_phi = new float[20];

	int partnum;
	int ijetnum=0;
	int nijet=0,ntaggedijet=0;

  //  TH2F *pt = new TH2F("jet_particle_pt","jet_particle_pt",1300,10,140,1300,10,140);
  //  TH2F *eta = new TH2F("jet_particle_eta","jet_particle_eta",60,-3,3,80,-4,4);
  //  TH2F *phi = new TH2F("jet_particle_phi","jet_particle_phi",60,-pi,pi,60,-pi,pi);
    TH1F *dr1D = new TH1F("dr","dr",70,0.,0.7);

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
				if ( dr > jet_R_size ) continue;
              //  if ( dr > 0.7) continue; //for _07finding
				if( p_pt[part]>leading_pt) {
					leading_pt = p_pt[part];
					leading_index = part; 
				}
			}//part
            nijet++;

            tagged_id = fabs( p_id[leading_index] );
            if( (mode==1 and (tagged_id==1 or tagged_id==2 or tagged_id==3)) or (mode==0 and tagged_id==21)){
                ntaggedijet++;
               
              //  pt->Fill(ijet_pt[ijet], p_pt[leading_index]);
              //  eta->Fill(ijet_eta[ijet], p_eta[leading_index]);
              //  phi->Fill(ijet_phi[ijet], p_phi[leading_index]);
                if(p_pt[leading_index]<ptcut[1] or ijet_pt[ijet]<ptcut[0]){
                    dr1D->Fill( getr(ijet_eta[ijet],p_eta[leading_index],ijet_phi[ijet],p_phi[leading_index]) );}
            }
            
            
		}//ijet
    }//ien
    cout<<"tagged: "<<ntaggedijet<<" /all jets in condition: "<<nijet<<endl;

	TString outpath="./cutteddR/jetparticle_cutteddR2_";

    if(mode==1) outpath.Append("gg2qqbar_");
    else outpath.Append("gg2gg_");

    if(inputpath.Contains("MIF")) outpath.Append("MIF_");
    else if(inputpath.Contains("MPI")) outpath.Append("MPI_");
    else if(inputpath.Contains("ISR")) outpath.Append("ISR_");
    else if(inputpath.Contains("FSR")) outpath.Append("FSR_");

    TString Rsize;
	Rsize.Form("R%02d",Rindex+2);
    outpath.Append(Rsize);
    outpath.Append(".root");
    TFile *outfile = new TFile(outpath, "RECREATE");
    
    //pt->Write();
    //eta->Write();
    //phi->Write();
    dr1D->Write();

    outfile->Write();
    outfile->Close();
return;
}
