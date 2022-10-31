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
//char* outputpath="./jetparticlegg2ggetaphi.root";

void partonmatch() //RMAX-1 -> R=0.7
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
	int ijetnum[RMAX]={0,0,0,0,0,0};

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
    
    TGraph *etaphi[15];
    for(int g=0; g<15; g++){ etaphi[g] = new TGraph(4);}

//R=0.7(Rindex=5) eta=-0.2~0.2
	float leading_pt, leading_pt_R02; 
    int leading_index, leading_index_R02;
    int tagged_id;
    float dr;
	int nijet=0, nmatched=0,nsameparton=0;
    int graphindex=0;

    const int entrysum=T->GetEntries();

	for (int ien=0; ien<entrysum; ien++){
		T->GetEntry(ien);

		for (int ijet=0; ijet<ijetnum[5]; ijet++){
            if(ijet_pt[5][ijet]>140 or ijet_pt[5][ijet]<10) continue;
            if(fabs(ijet_eta[5][ijet])>3) continue;
            nijet++;

            leading_pt = 0;
			for (int part=0; part<partnum; part++){
                
				dr = getr(ijet_eta[5][ijet],p_eta[part],ijet_phi[5][ijet],p_phi[part]);
				if ( dr > 0.7 ) continue;
				if( p_pt[part]>leading_pt) {
					leading_pt = p_pt[part];
					leading_index = part; 
				}
			}//part
             

            for (int sjet=0; sjet<ijetnum[0]; sjet++){
                if(ijet_pt[0][sjet]>140 or ijet_pt[0][sjet]<10) continue;
                if(fabs(ijet_eta[0][sjet])>3) continue;
                float rcone = getr(ijet_eta[0][sjet],ijet_eta[5][ijet],ijet_phi[0][sjet],ijet_phi[5][ijet]) ;
                if(rcone>0.2) continue;
                nmatched++;

                leading_pt_R02 = 0;
                for (int part=0; part<partnum; part++){
                
				    dr = getr(ijet_eta[0][sjet],p_eta[part],ijet_phi[0][sjet],p_phi[part]);
				    if ( dr > 0.2 ) continue;
                    

				    if( p_pt[part]>leading_pt_R02) {
					    leading_pt_R02 = p_pt[part];
					    leading_index_R02 = part; 
			    	}

                }  
                if( leading_index_R02==leading_index) nsameparton++;
                 else{

			        cout  
                      << "dR R02,R07=" << rcone
                      << ", entry #" <<ien
                      << ", jet R07#" <<ijet
                      << ", jet R02#" <<sjet
                      << endl
                      << "R07 leading mom id=" << p_id[leading_index]
                      << ", status="<< p_status[leading_index]
 				      << ", pT=" << p_pt[leading_index]
                      << ", dR mom,R07=" << getr(ijet_eta[5][ijet],p_eta[leading_index],ijet_phi[5][ijet],p_phi[leading_index])
                      << endl;
		        	cout 
                        << "R02 leading mom id=" << p_id[leading_index_R02]
                        << ", status="<< p_status[leading_index_R02]
 			        	<< ", pT=" << p_pt[leading_index_R02]
                      
                        << ", dR mom,R07=" <<  getr(ijet_eta[5][ijet],p_eta[leading_index_R02],ijet_phi[5][ijet],p_phi[leading_index_R02])
                        << endl << endl;
                    /*
                    etaphi[graphindex]->SetPoint(0,ijet_eta[5][ijet],ijet_phi[5][ijet]);
                    etaphi[graphindex]->SetPoint(1,p_eta[leading_index],p_phi[leading_index]);
                    etaphi[graphindex]->SetPoint(2,ijet_eta[0][sjet],ijet_phi[0][sjet]);
                    etaphi[graphindex]->SetPoint(3,p_eta[leading_index_R02],p_phi[leading_index_R02]);
                    graphindex++;
                    */
                     }
            } //sjet  


        

		}//ijet
    
	}//ien
cout<<
    nijet<<endl<<
    nmatched<<endl<<
    nsameparton<<endl;
/*
TFile *out = new TFile(outputpath, "recreate");
for(int g=0; g<graphindex; g++){ etaphi[g]->Write();}



out->Write();
out->Close();
*/
return;
}
