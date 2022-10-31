#include <iostream>
#include <vector>
#define pi 3.141592
#define RMAX 6
#define MOTHER 30
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
char* inputpath="./gg2qqbarfull.root";
void jetparticle2D(int Rindex=0, char* outpath="./jetparticlegg2qqR02_07finding.root") //RMAX-1 -> R=0.7
{
	float *p_pt=new float[500];
	float *p_eta=new float[500];
	float *p_phi=new float[500];
	bool *p_chg=new bool[500];
	int *p_id=new int[500];
    int *p_status=new int[500];

	float chgjet_pt[RMAX][10];
	float chgjet_eta[RMAX][10];
	float chgjet_phi[RMAX][10];
	float ijet_pt[RMAX][10];
	float ijet_eta[RMAX][10];
	float ijet_phi[RMAX][10];

	int partnum;
	int chgjetnum[RMAX]={0,0,0,0,0,0};
	int nchgjet=0,ntaggedchgjet=0;
	int ijetnum[RMAX]={0,0,0,0,0,0};
	int nijet=0,ntaggedijet=0;

    TH2F *pt = new TH2F("jet_particle_pt","jet_particle_pt",1300,10,140,1300,10,140);
    TH2F *eta = new TH2F("jet_particle_eta","jet_particle_eta",60,-3,3,80,-4,4);
    TH2F *phi = new TH2F("jet_particle_phi","jet_particle_phi",60,-pi,pi,60,-pi,pi);
    TH2F *chgpt = new TH2F("chgjet_particle_pt","chgjet_particle_pt",1300,10,140,1300,10,140);
    TH2F *chgeta = new TH2F("chgjet_particle_eta","chgjet_particle_eta",60,-3,3,80,-4,4);
    TH2F *chgphi = new TH2F("chgjet_particle_phi","chgjet_particle_phi",60,-pi,pi,60,-pi,pi);

    TH1F *dr1D = new TH1F("dr","dr",7,0.,0.7);

	TFile *f = new TFile(inputpath, "read");
	TTree *T = (TTree*)f ->Get("T");
	T->SetBranchAddress("p_pt",p_pt);
	T->SetBranchAddress("p_id",p_id);
    T->SetBranchAddress("p_status",p_status);
	T->SetBranchAddress("np",&partnum);
	T->SetBranchAddress("p_eta",p_eta);
	T->SetBranchAddress("p_phi",p_phi);
	T->SetBranchAddress("p_chg",p_chg);

	TString schgjet;
    TString sijet;

	for(int R=0;R<RMAX;R++){
		schgjet.Form("chgjetR%02d",R+2);
		T->SetBranchAddress(schgjet+"_pt",chgjet_pt[R]);
		T->SetBranchAddress(schgjet+"_eta",chgjet_eta[R]);
		T->SetBranchAddress(schgjet+"_phi",chgjet_phi[R]);
		T->SetBranchAddress("n"+schgjet,&chgjetnum[R]);

		sijet.Form("ijetR%02d",R+2);
		T->SetBranchAddress(sijet+"_pt",ijet_pt[R]);
		T->SetBranchAddress(sijet+"_eta",ijet_eta[R]);
		T->SetBranchAddress(sijet+"_phi",ijet_phi[R]);
		T->SetBranchAddress("n"+sijet,&ijetnum[R]);

	}
    
    
	float jet_R_size = 0.1*(2.0+(float)Rindex); //R=0.7(Rindex=5) eta=-0.2~0.2
	float leading_pt; 
    int leading_index;
    int tagged_id;
    float dr;

	const int entrysum=T->GetEntries();

	for (int ien=0; ien<entrysum; ien++){
		T->GetEntry(ien);

		//if ( ijetnum[Rindex]<1 ) continue;
        
       // if( ijetnum[Rindex] < chgjetnum[Rindex] ) cout<<ien<<endl;
        
		for (int ijet=0; ijet<ijetnum[Rindex]; ijet++){
            if(ijet_pt[Rindex][ijet]>140 or ijet_pt[Rindex][ijet]<10) continue;
            if(fabs(ijet_eta[Rindex][ijet])>3) continue;
            leading_pt=0;
			for (int part=0; part<partnum; part++){
                
				dr = getr(ijet_eta[Rindex][ijet],p_eta[part],ijet_phi[Rindex][ijet],p_phi[part]);
			//	if ( dr > jet_R_size ) continue;
                if ( dr > 0.7) continue; //for _07finding
				if( p_pt[part]>leading_pt) {
					leading_pt = p_pt[part];
					leading_index = part; 
				}
			}//part
            nijet++;

            tagged_id = fabs( p_id[leading_index] );
            if( tagged_id==1 or tagged_id==2 or tagged_id==3){
            //if( tagged_id == 21){
                ntaggedijet++;
               
                pt->Fill(ijet_pt[Rindex][ijet], p_pt[leading_index]);
                eta->Fill(ijet_eta[Rindex][ijet], p_eta[leading_index]);
                phi->Fill(ijet_phi[Rindex][ijet], p_phi[leading_index]);

                dr1D->Fill( getr(ijet_eta[Rindex][ijet],p_eta[leading_index],ijet_phi[Rindex][ijet],p_phi[leading_index]) );
            }
            else{
            
           
			cout 
				<< "jet pT=" << ijet_pt[Rindex][ijet] 
				<< ", eta=" << ijet_eta[Rindex][ijet] 
				<< endl;
			cout 
                << "leading mom id=" << p_id[leading_index]
                << ", status="<< p_status[leading_index]
 				<< ", pT=" << p_pt[leading_index]
				<< ", eta=" << p_eta[leading_index]
                << endl;
            }
            
		}//ijet
       
        for (int chgjet=0; chgjet<chgjetnum[Rindex]; chgjet++){
            if(chgjet_pt[Rindex][chgjet]>140 or chgjet_pt[Rindex][chgjet]<10) continue;
            if(fabs(chgjet_eta[Rindex][chgjet])>3) continue;
            leading_pt=0;
			for (int part=0; part<partnum; part++){
                
				dr = getr(chgjet_eta[Rindex][chgjet],p_eta[part],chgjet_phi[Rindex][chgjet],p_phi[part]);
			//	if ( dr > jet_R_size ) continue;
                if ( dr > 0.7) continue; //for _07finding
				if( p_pt[part]>leading_pt) {
					leading_pt = p_pt[part];
					leading_index = part; 
				}
			}//part
            nchgjet++;

            tagged_id = fabs( p_id[leading_index] );
            if( tagged_id==1 or tagged_id==2 or tagged_id==3){
           // if( tagged_id == 21){
                ntaggedchgjet++;
                chgpt->Fill(chgjet_pt[Rindex][chgjet], p_pt[leading_index]);
                chgeta->Fill(chgjet_eta[Rindex][chgjet], p_eta[leading_index]);
                chgphi->Fill(chgjet_phi[Rindex][chgjet], p_phi[leading_index]);
            }
         /*   else{
            
			cout 
				<< "jet pT=" << chgjet_pt[Rindex][chgjet] 
				<< ", eta=" << chgjet_eta[Rindex][chgjet] 
				<< endl;
			cout 
                << "leading mom id=" << p_id[leading_index]
                << ", status="<< p_status[leading_index]
 				<< ", pT=" << p_pt[leading_index]
				<< ", eta=" << p_eta[leading_index]
                << endl;             }
          */
		}//chgjet

    
	}//ien
    cout<<"tagged: "<<ntaggedijet<<" /all jets in condition: "<<nijet<<endl;
    cout<<"tagged(chg): "<<ntaggedchgjet<<" /all jets in condition(chg): "<<nchgjet<<endl;
    TFile *outfile = new TFile(outpath, "RECREATE");
    pt->Write();
    eta->Write();
    phi->Write();

    dr1D->Write();

    chgpt->Write();
    chgeta->Write();
    chgphi->Write();

    outfile->Write();
    outfile->Close();
return;
}
