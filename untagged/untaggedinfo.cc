#include <iostream>
#include <fstream>
#include <vector>
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
void untaggedinfo(int Rindex=0, TString inputpath="../gg2qqbarfull_MIF.root") //RMAX-1 -> R=0.7
{
    int mode=inputpath.Contains("qq"); //if quark jet mode =1 

    float *p_pt=new float[2000];
	float *p_eta=new float[2000];
	float *p_phi=new float[2000];
	bool *p_chg=new bool[2000];
	int *p_id=new int[2000];
    int *p_status=new int[2000];

	float *ijet_pt=new float[20];
	float *ijet_eta=new float[20];
	float *ijet_phi=new float[20];

	int partnum;
	int ijetnum=0;
	int nijet=0,ntaggedijet=0;

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
    int qgconf=0;
    
    
    TString outpath="./untagged_";
    if(mode==1) outpath.Append("gg2qqbar");
    else outpath.Append("gg2gg");

    if(inputpath.Contains("MIF")) outpath.Append("MIF");

    TString Rsize;
	Rsize.Form("R%02d",Rindex+2);
    outpath.Append(Rsize);
    outpath.Append(".txt");

   // ofstream fout(outpath);
    ofstream fout("test.txt");

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
            if( (mode==1 and (tagged_id==1 or tagged_id==2 or tagged_id==3)) or (mode==0 and tagged_id==21))  ntaggedijet++;
            else if( (mode==0 and (tagged_id==1 or tagged_id==2 or tagged_id==3)) or (mode==1 and tagged_id==21) ){
                qgconf++;
            }

            else{
          
           
			fout 
				<< "jet pT=" << ijet_pt[ijet] 
				<< ", eta=" << ijet_eta[ijet] 
				<< endl;
			fout 
                << "leading mom id=" << p_id[leading_index]
                << ", status="<< p_status[leading_index]
 				<< ", pT=" << p_pt[leading_index]
				<< ", eta=" << p_eta[leading_index]
                << endl;

                
            }
            
		}//ijet
    }//ien
    
    cout
        <<endl
        <<"numbers:"<<endl
        <<"untagged: "<<nijet-ntaggedijet<<" /all jets in condition: "<<nijet<<endl
        <<"confused with light quark or gluon: "<<qgconf<<endl;
    fout
        <<endl
        <<"numbers:"<<endl
        <<"untagged: "<<nijet-ntaggedijet<<" /all jets in condition: "<<nijet<<endl
        <<"confused with light quark or gluon: "<<qgconf<<endl;

    fout.close();
    
    return;
}
