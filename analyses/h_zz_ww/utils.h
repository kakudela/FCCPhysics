

namespace FCCAnalyses {

bool is_ww_hadronic(Vec_mc mc, Vec_i ind) {
   int l1 = 0;
   int l2 = 0;
   //cout << "*********" << endl;
   for(size_t i = 0; i < mc.size(); ++i) {
        auto & p = mc[i];
        if(std::abs(p.PDG) == 24) {
            int ds = p.daughters_begin;
            int de = p.daughters_end;
            for(int k=ds; k<de; k++) {
                int pdg = abs(mc[ind[k]].PDG);
                if(pdg == 24) continue;
                //std::cout << "W " << pdg << endl;
                if(l1 == 0) l1 = pdg;
                else l2 = pdg;
            }
        }
   }
   if((l1 < 6 && l2 < 6)) {
       //std::cout << "HADRONIC-----------" << l1 << " " << l2 << endl;
       return true;
   }
   return false;
}

bool is_zz_hadronic(Vec_mc mc, Vec_i ind) {
   int l1 = 0;
   int l2 = 0;
   //cout << "*********" << endl;
   for(size_t i = 0; i < mc.size(); ++i) {
        auto & p = mc[i];
        if(std::abs(p.PDG) == 23) {
            int ds = p.daughters_begin;
            int de = p.daughters_end;
            for(int k=ds; k<de; k++) {
                int pdg = abs(mc[ind[k]].PDG);
                if(pdg == 24) continue;
                //std::cout << "W " << pdg << endl;
                if(l1 == 0) l1 = pdg;
                else l2 = pdg;
            }
        }
   }
   if((l1 < 6 && l2 < 6)) {
       //std::cout << "HADRONIC-----------" << l1 << " " << l2 << endl;
       return true;
   }
   return false;
}

Vec_tlv pairing_WW(Vec_tlv J, float target=80.385) {

    int nJets6=6;
    float d1W[6][6], dW1[6][6], d2W[6][6],dW2[6][6], d1Z[6][6],d2Z[6][6], d1H[6][6][6][6],d2H[6][6][6][6];

    for(int i=0;i<6;i++) {
        for(int j=0;j<6;j++) {
            d1W[i][j]=-1000.;
            d2W[i][j]=-1000.;
            d1Z[i][j]=-1000.;
            d2Z[i][j]=-1000.;
            dW1[i][j]=-1000.;
            dW2[i][j]=-1000.;
        }
    }

    for(int i=0;i<6;i++){
        for(int j=0;j<6;j++){
            for(int k=0;k<6;k++){
                for(int s=0;s<6;s++) {
                    d1H[i][j][k][s]=-1000.;
                    d2H[i][j][k][s]=-1000.;
                }
            }
        }
    }

    double DWmin_1=1000.,DWmin_2=1000.,chi2_1=10000000000., chi2_2=10000000000.,chi2M=10000000000.,chi2min=100000000.,chi2First=10000000.,chi2Second=100000000;float dmin=1000000.,DZmin=1000000.;

    int iW1=-1.,jW1=-1.,iW2=-1.,jW2=-1.,iZ=-1.,jZ=-1.;
    int iW1_Z=-1.,jW1_Z=-1.,iW2_Z=-1.,jW2_Z=-1.,iZ_Z=-1.,jZ_Z=-1.;
    int iW1_1=-1.,jW1_1=-1.,iW2_1=-1.,jW2_1=-1.,iZ_1=-1.,jZ_1=-1.;
    int iW1_2=-1.,jW1_2=-1.,iW2_2=-1.,jW2_2=-1.,iZ_2=-1.,jZ_2=-1.;
    int iW1_3=-1.,jW1_3=-1.,iW2_3=-1.,jW2_3=-1.,iZ_3=-1.,jZ_3=-1.;


    for(int i=0;i<nJets6;i++) {
        for(int j=i+1;j<nJets6;j++) {
            for(int k=i+1;k<nJets6;k++) {
                if(!(k==i)&&!(k==j)) {
                    for(int s=k+1;s<nJets6;s++) {
                        if(!(s==i)&&!(s==j)&&!(s==k)) {
                            for(int l=0;l<nJets6;l++) {
                                if(!(l==i)&&!(l==j)&&!(l==k)&&!(l==s)) {
                                    for(int m=l+1;m<nJets6;m++) {
                                        if(!(m==i)&&!(m==j)&&!(m==k)&&!(m==s)&&!(m==l)) {

                                            d1W[i][j]      =fabs( target -((J[i]+J[j]).M()));
                                            dW1[l][m]      =fabs( target -((J[l]+J[m]).M()));
                                            d1Z[k][s]      =fabs( 91.19  -((J[k]+J[s]).M()));
                                            d1H[i][j][l][m]=fabs(125.0   -((J[i]+J[j]+J[l]+J[m]).M()));
                                      
                                            chi2First=(d1W[i][j])*(d1W[i][j])+(d1Z[k][s])*(d1Z[k][s]); 
                                            chi2_1   =(d1W[i][j])*(d1W[i][j])+(d1Z[k][s])*(d1Z[k][s])+d1H[i][j][l][m]*d1H[i][j][l][m]; 

                                            if(d1Z[k][s]<DZmin) {
                                                DZmin=d1Z[k][s];
                                                iZ_Z=k;
                                                jZ_Z=s;
                                                if(d1W[i][j]<dW1[l][m]) { iW1_Z=i; jW1_Z=j; iW2_Z=l;jW2_Z=m;}
                                                else                    { iW1_Z=l; jW1_Z=m; iW2_Z=i;jW2_Z=j;}
                                            }
                                            if(chi2First<chi2min) {chi2min=chi2First;iW1_1=i;jW1_1=j; iZ_1=k; jZ_1=s;iW2_1=l;jW2_1=m;}
                                            if(chi2_1<chi2M)      {chi2M=chi2_1;     iW1_2=i;jW1_2=j; iZ_2=k; jZ_2=s;iW2_2=l;jW2_2=m;}
                                            if( DWmin_1<dmin)     {dmin= DWmin_1;    iW1_3=i;jW1_3=j; iZ_3=k; jZ_3=s;iW2_3=l;jW2_3=m;}

                                            d2W[k][s]      =fabs( target -((J[k]+J[s]).M()));
                                            dW2[l][m]      =fabs( target -((J[l]+J[m]).M()));
                                            d2Z[i][j]      =fabs( 91.19  -((J[i]+J[j]).M()));
                                            d2H[k][s][l][m]=fabs(125.0   -((J[k]+J[s]+J[l]+J[m]).M()));
                                      
                                          
                                            chi2Second=(d2W[k][s])*(d2W[k][s])+(d2Z[i][j])*(d2Z[i][j]); 
                                            chi2_2=(d2W[k][s])*(d2W[k][s])+(d2Z[i][j])*(d2Z[i][j])+(d2H[k][s][l][m]*d2H[k][s][l][m]); 

                                            if(d2Z[i][j]<DZmin) {
                                                DZmin=d2Z[i][j];
                                                iZ_Z=i;
                                                jZ_Z=j; 

                                                if(d2W[k][s]<dW2[l][m]) {iW1_Z=k; jW1_Z=s; iW2_Z=l; jW2_Z=m;}   
                                                else {iW1_Z=l; jW1_Z=m; iW2_Z=k; jW2_Z=s;}
                                            }

                                            if(chi2Second<chi2min) {chi2min=chi2Second;iW1_1=k;jW1_1=s;iZ_1=i; jZ_1=j;iW2_1=l;jW2_1=m;}
                                            if(chi2_2<chi2M)       {chi2M  =chi2_2;    iW1_2=k;jW1_2=s;iZ_2=i; jZ_2=j;iW2_2=l;jW2_2=m;}
                                            if( DWmin_2<dmin)      {dmin   =DWmin_2;   iW1_3=k;jW1_3=s; iZ_3=i;jZ_3=j;iW2_3=l;jW2_3=m;}
                                        }
                                    } // m
                                }
                            }// l
                        }
                    } // s
                }
            }// k
        } // j
    
    }// i

    TLorentzVector JHiggs,JW1,JW2,JZ;
    JW1=J[iW1_2]+J[jW1_2];
    JW2=J[iW2_2]+J[jW2_2];
    JZ=J[iZ_2]+J[jZ_2];
    JHiggs=JW1+JW2;

    Vec_tlv out;
    out.push_back(JW1);
    out.push_back(JW2);
    out.push_back(JZ);
    out.push_back(JHiggs);
    return out;

}




}