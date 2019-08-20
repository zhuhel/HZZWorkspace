#include "HZZWorkspace/IntegralHelper.h"

#include <TMath.h>
#include <math.h>

using namespace std;

namespace IntegralHelper{
    double L_gg (double x)
    {
        double m[250]={100, 101.2103718, 102.4353935, 103.6752426, 104.9300984, 106.2001427, 107.4855592, 108.7865341, 110.1032556, 111.4359143, 112.7847031, 114.1498173, 115.5314544, 116.9298145, 118.3451, 119.7775156, 121.2272689, 122.6945695, 124.1796299, 125.6826651, 127.2038925, 128.7435325, 130.3018079, 131.8789442, 133.4751696, 135.0907154, 136.7258153, 138.3807059, 140.0556269, 141.7508206, 143.4665325, 145.2030109, 146.9605071, 148.7392756, 150.5395738, 152.3616623, 154.2058048, 156.0722683, 157.961323, 159.8732422, 161.8083028, 163.7667848, 165.7489717, 167.7551504, 169.7856114, 171.8406484, 173.9205591, 176.0256444, 178.1562091, 180.3125615, 182.4950139, 184.703882, 186.9394856, 189.2021483, 191.4921977, 193.8099652, 196.1557862, 198.5300005, 200.9329515, 203.3649872, 205.8264596, 208.3177249, 210.8391438, 213.3910812, 215.9739066, 218.5879938, 221.2337211, 223.9114716, 226.6216328, 229.364597, 232.1407613, 234.9505275, 237.7943023, 240.6724974, 243.5855293, 246.5338198, 249.5177955, 252.5378884, 255.5945357, 258.6881798, 261.8192684, 264.9882549, 268.1955979, 271.4417617, 274.7272161, 278.0524367, 281.4179049, 284.8241077, 288.2715382, 291.7606955, 295.2920846, 298.8662165, 302.4836088, 306.144785, 309.850275, 313.6006152, 317.3963485, 321.2380242, 325.1261985, 329.0614342, 333.0443009, 337.075375, 341.1552401, 345.2844868, 349.4637127, 353.6935228, 357.9745293, 362.3073519, 366.6926177, 371.1309616, 375.6230259, 380.1694609, 384.7709247, 389.4280833, 394.1416108, 398.9121895, 403.74051, 408.6272711, 413.5731802, 418.5789531, 423.6453146, 428.7729978, 433.962745, 439.2153075, 444.5314456, 449.9119286, 455.3575355, 460.8690545, 466.4472834, 472.0930295, 477.8071102, 483.5903525, 489.4435936, 495.3676806, 501.3634711, 507.4318329, 513.5736445, 519.7897948, 526.0811837, 532.4487217, 538.8933307, 545.4159433, 552.0175038, 558.6989678, 565.4613023, 572.3054862, 579.2325101, 586.2433768, 593.3391011, 600.52071, 607.789243, 615.1457523, 622.5913028, 630.1269721, 637.7538509, 645.4730434, 653.2856668, 661.192852, 669.1957435, 677.2954998, 685.4932932, 693.7903104, 702.1877524, 710.6868346, 719.2887873, 727.9948556, 736.8062997, 745.7243951, 754.7504325, 763.8857186, 773.1315755, 782.4893418, 791.9603717, 801.5460364, 811.2477232, 821.0668365, 831.0047976, 841.0630449, 851.2430345, 861.5462397, 871.9741521, 882.5282809, 893.2101539, 904.0213173, 914.963336, 926.0377938, 937.2462937, 948.5904581, 960.0719291, 971.6923685, 983.4534585, 995.3569014, 1007.40442, 1019.597759, 1031.938682, 1044.428976, 1057.07045, 1069.864932, 1082.814275, 1095.920353, 1109.185063, 1122.610326, 1136.198084, 1149.950305, 1163.868979, 1177.95612, 1192.213768, 1206.643987, 1221.248865, 1236.030516, 1250.991081, 1266.132723, 1281.457636, 1296.968037, 1312.666172, 1328.554313, 1344.634759, 1360.909838, 1377.381906, 1394.053348, 1410.926576, 1428.004033, 1445.28819, 1462.78155, 1480.486645, 1498.406037, 1516.542321, 1534.89812, 1553.476094, 1572.27893, 1591.30935, 1610.570109, 1630.063994, 1649.793828, 1669.762467, 1689.9728, 1710.427754, 1731.130288, 1752.0834, 1773.290123, 1794.753525, 1816.476715, 1838.462836, 1860.715071, 1883.236641, 1906.030805, 1929.100864, 1952.450156, 1976.082061, 2000};

        double y[250]={-3.23043785351128, -3.27176435438805, -3.31311508683479, -3.35458693019852, -3.39610963979821, -3.4376976364921, -3.47943707392032, -3.52119221956324, -3.56307897700638, -3.60499763983444, -3.64701871121563, -3.68912096047425, -3.73129669953724, -3.77357191793491, -3.81595251927541, -3.85834223360931, -3.90082509174404, -3.94340953397927, -3.98609611826638, -4.02888234994357, -4.06523525030659, -4.11460185225372, -4.15760120132163, -4.20069018093465, -4.24384461380463, -4.28710670991824, -4.33049077366799, -4.37385690518053, -4.41734750068601, -4.46104137575284, -4.50463411677953, -4.54838561101719, -4.59222361468643, -4.63615101392609, -4.68015534096959, -4.7242591587092, -4.76844596947466, -4.81271369763732, -4.85702872236897, -4.90151990916171, -4.9460627536936, -4.99067949485503, -5.03538439382575, -5.08017655593957, -5.12492171505367, -5.17003227283999, -5.21501043361317, -5.26023636377181, -5.30549371202873, -5.35083672677354, -5.39626768587563, -5.44175362520327, -5.48731907241776, -5.53304193951674, -5.57876883546556, -5.62472418234513, -5.67068043786739, -5.71673096154985, -5.76290091325503, -5.80912273418585, -5.85546440256077, -5.90189317224423, -5.94831324960149, -5.99507379581708, -6.04184782999741, -6.08859655707755, -6.13549045248636, -6.18248650768751, -6.22963321883022, -6.27680843258403, -6.32414953473881, -6.37150288077392, -6.41899659423632, -6.46668198723454, -6.51427413592681, -6.5620984728035, -6.60997883902126, -6.65795716116356, -6.70605284149434, -6.75430670128519, -6.80250956141336, -6.85094546784149, -6.89944473144302, -6.94805327081578, -6.99676563299571, -7.04558949388343, -7.0945133363071, -7.14353731263965, -7.19265687600754, -7.24191387758909, -7.29124225923596, -7.34070031723859, -7.39021606720835, -7.43981965984523, -7.48972026021728, -7.53953081479324, -7.58958701395148, -7.63974491817862, -7.68994158557135, -7.74027006301905, -7.79072341490147, -7.84150787440376, -7.89194893458788, -7.94276155084723, -7.99366969208403, -8.04471514170044, -8.09585165410403, -8.14704976219758, -8.19846055810363, -8.24997069842896, -8.30158256628179, -8.35333753303873, -8.40521020669396, -8.45728733043473, -8.50925386146863, -8.56149779108879, -8.61390927529071, -8.66631059039709, -8.71867621250331, -8.77163101979388, -8.82448521084824, -8.87746739205804, -8.93049101959129, -8.98381108673312, -9.03714433974831, -9.09068722118679, -9.14433994993664, -9.19807704094381, -9.25198312426763, -9.30608120857984, -9.36018545747532, -9.41450055453408, -9.46894785117661, -9.52356059969548, -9.57826077892407, -9.63310741809194, -9.68811590856508, -9.74320596615493, -9.79854938995681, -9.85398102920368, -9.90959427023927, -9.9652981394187, -10.0211577225121, -10.077182698512, -10.1333670852194, -10.1896775371049, -10.2461485535147, -10.3027662795212, -10.3595827119735, -10.4164864524903, -10.4736059151353, -10.5308404633659, -10.5882445659924, -10.6458123588791, -10.7035021021541, -10.7614392561033, -10.8194951469975, -10.8777286273315, -10.9360673092818, -10.9946796919289, -11.0533674555005, -11.1124401896728, -11.1713463499646, -11.2305880296806, -11.2899885133963, -11.3496006932505, -11.4093507353796, -11.4693025346555, -11.5294286051537, -11.5897610267552, -11.650230887667, -11.710903802217, -11.7717769068998, -11.8328398960792, -11.8940854046568, -11.9555102381239, -12.0171156025155, -12.0789438782115, -12.1409688740517, -12.2028528388755, -12.265663538391, -12.3282157150603, -12.3910326449146, -12.4540635034805, -12.5173288537475, -12.5807287356695, -12.6443152656971, -12.7082713537048, -12.7723587223114, -12.8368500921973, -12.9011813830857, -12.9659162005461, -13.0307540532351, -13.0960810651539, -13.161496417241, -13.2271497068275, -13.2930310639056, -13.3591517618672, -13.4255140264532, -13.4921038317478, -13.5589232944335, -13.6260814133557, -13.6932397045324, -13.7609317830336, -13.8287386277735, -13.8968358813127, -13.965171589197, -14.0338237966695, -14.1026426783341, -14.1716452681413, -14.2409644438255, -14.3108210120341, -14.3807614637448, -14.4509830242569, -14.5214603272422, -14.5922903654874, -14.6632598783109, -14.7346896180453, -14.8063455552311, -14.8782967063161, -14.9505357072477, -15.0230896809139, -15.0959623150496, -15.1691151589447, -15.2425963115727, -15.3163955876211, -15.3905316389688, -15.464966827505, -15.5397410977859, -15.614856420441, -15.6902787106767, -15.7660797694148, -15.8422025453579, -15.9186890237303, -15.9955246629723, -16.0727326432627, -16.1502740813364, -16.2282012866903, -16.306499698048, -16.3851738515752, -16.4642289117701, -16.5436660657705, -16.6235233051784, -16.7037632853722, -16.7844029591119, -16.8654492265505, -16.946940503461, -17.0287860169476, -17.1111246487548, -17.1938575623246}; //CT10nnlo


        for (int i=0; i<248; i++) {if  (x>=m[i+0] && x<=m[i+1]) return exp( y[i+1]-((y[i]-y[i+1]-(y[i+2]-y[i+1])*((m[i]*m[i]-m[i+1]*m[i+1])/(m[i+2]*m[i+2]-m[i+1]*m[i+1])))/(m[i]-m[i+1]-(m[i+2]-m[i+1])*((m[i]*m[i]-m[i+1]*m[i+1])/(m[i+2]*m[i+2]-m[i+1]*m[i+1]))))*m[i+1] - ((y[i]-y[i+1]-(y[i+2]-y[i+1])*((m[i]-m[i+1])/(m[i+2]-m[i+1])))/(m[i]*m[i]-m[i+1]*m[i+1]-(m[i+2]*m[i+2]-m[i+1]*m[i+1])*((m[i]-m[i+1])/(m[i+2]-m[i+1]))))*m[i+1]*m[i+1] + ((y[i]-y[i+1]-(y[i+2]-y[i+1])*((m[i]*m[i]-m[i+1]*m[i+1])/(m[i+2]*m[i+2]-m[i+1]*m[i+1])))/(m[i]-m[i+1]-(m[i+2]-m[i+1])*((m[i]*m[i]-m[i+1]*m[i+1])/(m[i+2]*m[i+2]-m[i+1]*m[i+1]))))*x  + ((y[i]-y[i+1]-(y[i+2]-y[i+1])*((m[i]-m[i+1])/(m[i+2]-m[i+1])))/(m[i]*m[i]-m[i+1]*m[i+1]-(m[i+2]*m[i+2]-m[i+1]*m[i+1])*((m[i]-m[i+1])/(m[i+2]-m[i+1]))))*x*x );}
        if  (x>=m[248] && x<=m[249]) {
            int i = 247;
            return exp( y[i+1]-((y[i]-y[i+1]-(y[i+2]-y[i+1])*((m[i]*m[i]-m[i+1]*m[i+1])/(m[i+2]*m[i+2]-m[i+1]*m[i+1])))/(m[i]-m[i+1]-(m[i+2]-m[i+1])*((m[i]*m[i]-m[i+1]*m[i+1])/(m[i+2]*m[i+2]-m[i+1]*m[i+1]))))*m[i+1] - ((y[i]-y[i+1]-(y[i+2]-y[i+1])*((m[i]-m[i+1])/(m[i+2]-m[i+1])))/(m[i]*m[i]-m[i+1]*m[i+1]-(m[i+2]*m[i+2]-m[i+1]*m[i+1])*((m[i]-m[i+1])/(m[i+2]-m[i+1]))))*m[i+1]*m[i+1] + ((y[i]-y[i+1]-(y[i+2]-y[i+1])*((m[i]*m[i]-m[i+1]*m[i+1])/(m[i+2]*m[i+2]-m[i+1]*m[i+1])))/(m[i]-m[i+1]-(m[i+2]-m[i+1])*((m[i]*m[i]-m[i+1]*m[i+1])/(m[i+2]*m[i+2]-m[i+1]*m[i+1]))))*x  + ((y[i]-y[i+1]-(y[i+2]-y[i+1])*((m[i]-m[i+1])/(m[i+2]-m[i+1])))/(m[i]*m[i]-m[i+1]*m[i+1]-(m[i+2]*m[i+2]-m[i+1]*m[i+1])*((m[i]-m[i+1])/(m[i+2]-m[i+1]))))*x*x );}

        if  (x>m[249]) return exp(y[249]);
        else if  (x<m[0])  return exp(y[0]);
        else return 1;
    }

    Double_t getHiggsXS(Double_t mH, Double_t width)
    {
        Double_t Norm = 6.88734e-08;

        Double_t w_fr = width / mH;

        Double_t min_x = mH - 10*width;
        if (min_x<230) min_x=230;
        Double_t max_x = mH + 10*width;
        Double_t n_steps  = 20000;
        Double_t step  = (max_x - min_x) / n_steps;


        Double_t mH_cmplx = mH * sqrt(1+w_fr*w_fr);

        Double_t result = 0;
        for (int i_step=0; i_step<n_steps; i_step++)
        {
            Double_t x = min_x + step*(i_step+0.5);

            Double_t A = 1 + w_fr*w_fr;
            Double_t B = std::pow((x*x-mH_cmplx*mH_cmplx),2) + std::pow((x*x*w_fr),2);
            Double_t Propag = A/B;

            Double_t G_gg = std::pow(x,3);
            Double_t t_t = std::pow(x/(2*173.34),2);
            Double_t t_b = std::pow(x/(2*4.18),2);
            Double_t A_2 = 0;
            Double_t pi = TMath::Pi();
            if (t_t<=1 && t_b<=1) {
                A_2 = pow( 1/t_t + (t_t-1)/(t_t*t_t)*pow(TMath::ASin(sqrt(t_t)),2) + 1/t_b + (t_b-1)/(t_b*t_b)*pow(TMath::ASin(sqrt(t_b)),2) ,2);
            } else if (t_t<=1 && t_b>1) {
                Double_t log_b = TMath::Log((1+sqrt(1-1/t_b)) / (1-sqrt(1-1/t_b)));
                A_2 = pow( 1/t_t + (t_t-1)/(t_t*t_t)*pow(TMath::ASin(sqrt(t_t)),2) + 1/t_b - 0.25*(t_b-1)/(t_b*t_b)*(log_b*log_b - pi*pi) ,2) + pow( 0.5*(t_b-1)/(t_b*t_b)*pi*log_b ,2);
            } else if (t_t>1 && t_b>1) {
                Double_t log_b = TMath::Log((1+sqrt(1-1/t_b)) / (1-sqrt(1-1/t_b)));
                Double_t log_t = TMath::Log((1+sqrt(1-1/t_t)) / (1-sqrt(1-1/t_t)));
                A_2 = pow( 1/t_b - 0.25*(t_b-1)/(t_b*t_b)*(log_b*log_b - pi*pi) + 1/t_t - 0.25*(t_t-1)/(t_t*t_t)*(log_t*log_t - pi*pi) ,2) + 0.25*pow( (t_b-1)/(t_b*t_b)*pi*log_b + (t_t-1)/(t_t*t_t)*pi*log_t ,2);
            }
            G_gg *= A_2;

            Double_t G_zz = 1;
            G_zz *= std::pow(x,3) * std::pow((1-std::pow(2*91.2/x,2)),1/2) * (1-std::pow(2*91.2/x,2)+0.75*std::pow(2*91.2/x,4)) * (x>2*91.2);

            Double_t result_tmp =  Norm * L_gg(x) * std::pow(x,1) * G_gg * Propag * G_zz ;
            //	if (i_step%100==0) std::cout<<"i_step = "<<i_step<<", x = "<<x<<", result_tmp = "<<result_tmp<< endl;
            if (result_tmp>0) result+= result_tmp;
        }
        result*=step/2.;


        //    std::cout<<"[get_stat::getHiggsXS] Higgs cross section computed: mH = "<<mH<<" GeV, Width = "<<width<<" GeV, xs = "<<result<<" fb^-1.\n";
        return result;
    }

    Double_t getHiggsXSTF1(Double_t *x, Double_t *par)
    {
        Double_t mH = x[0];
        Double_t width = mH * par[0];
        //    std::cout<<"[get_stat::getHiggsXSTF1] : mH = "<<mH<<" GeV, Width = "<<width<<std::endl;
        return getHiggsXS(mH, width);
    }

    Double_t getSignalIntegral(Double_t mH, Double_t width, Double_t acc_0, Double_t acc_1, Double_t acc_2)
    {
        Double_t Norm = 6.88734e-08;

        Double_t w_fr = width / mH;

        Double_t min_x = mH - 6*width;
        if (min_x<240) min_x=240;
        Double_t max_x = mH + 6*width;
        Double_t n_steps  = 3000;
        Double_t step  = (max_x - min_x) / n_steps;


        Double_t mH_cmplx = mH * sqrt(1+w_fr*w_fr);

        Double_t result = 0;
        for (int i_step=0; i_step<n_steps; i_step++)
        {
            Double_t x = min_x + step*(i_step+0.5);

            Double_t A = 1 + w_fr*w_fr;
            Double_t B = std::pow((x*x-mH_cmplx*mH_cmplx),2) + std::pow((x*x*w_fr),2);
            Double_t Propag = A/B;

            Double_t G_gg = std::pow(x,3);
            Double_t t_t = std::pow(x/(2*173.34),2);
            Double_t t_b = std::pow(x/(2*4.18),2);
            Double_t A_2 = 0;
            Double_t pi = TMath::Pi();
            if (t_t<=1 && t_b<=1) {
                A_2 = pow( 1/t_t + (t_t-1)/(t_t*t_t)*pow(TMath::ASin(sqrt(t_t)),2) + 1/t_b + (t_b-1)/(t_b*t_b)*pow(TMath::ASin(sqrt(t_b)),2) ,2);
            } else if (t_t<=1 && t_b>1) {
                Double_t log_b = TMath::Log((1+sqrt(1-1/t_b)) / (1-sqrt(1-1/t_b)));
                A_2 = pow( 1/t_t + (t_t-1)/(t_t*t_t)*pow(TMath::ASin(sqrt(t_t)),2) + 1/t_b - 0.25*(t_b-1)/(t_b*t_b)*(log_b*log_b - pi*pi) ,2) + pow( 0.5*(t_b-1)/(t_b*t_b)*pi*log_b ,2);
            } else if (t_t>1 && t_b>1) {
                Double_t log_b = TMath::Log((1+sqrt(1-1/t_b)) / (1-sqrt(1-1/t_b)));
                Double_t log_t = TMath::Log((1+sqrt(1-1/t_t)) / (1-sqrt(1-1/t_t)));
                A_2 = pow( 1/t_b - 0.25*(t_b-1)/(t_b*t_b)*(log_b*log_b - pi*pi) + 1/t_t - 0.25*(t_t-1)/(t_t*t_t)*(log_t*log_t - pi*pi) ,2) + 0.25*pow( (t_b-1)/(t_b*t_b)*pi*log_b + (t_t-1)/(t_t*t_t)*pi*log_t ,2);
            }
            G_gg *= A_2;

            Double_t G_zz = 1;
            G_zz *= std::pow(x,3) * std::pow((1-std::pow(2*91.2/x,2)),1/2) * (1-std::pow(2*91.2/x,2)+0.75*std::pow(2*91.2/x,4)) * (x>2*91.2);

            Double_t Acc =  acc_0 + acc_1*x + acc_2*x*x ;

            Double_t result_tmp =  Norm * L_gg(x) * std::pow(x,1) * G_gg * Propag * G_zz * Acc ;

            if (result_tmp>0) result+= result_tmp;
        }
        result*=step;


        //std::cout<<"[get_stat::getSignalIntergal] Signal integral computed: m_H = "<<mH<<" GeV, Width = "<<width<<" GeV, acc = "<<acc_0<<" "<<acc_1<<" "<<acc_2<<" "<<", xs = "<<result<<"\n";
        return result;
    }

    Double_t getSignalIntegralTF1(Double_t *x, Double_t *par)
    {
        /*    size_t x_size = sizeof(x)/sizeof(x[0]);
              size_t par_size = sizeof(par)/sizeof(par[0]);
              if( x_size < 2) {
              cerr <<"getSignalIntegralTF2: Not enough dimentions, at least 2, but only"<< x_size <<"given"<<endl;
              return 0;
              }
              if(par_size < 3){
              cerr <<"getSignalIntegralTF2: Not enough parameters, at least 3, but only"<< par_size <<"given"<<endl;
              return 0;
              }
              */
        Double_t mH = x[0];
        Double_t width = mH * par[0];
        Double_t acc_0 = par[1];
        Double_t acc_1 = par[2];
        Double_t acc_2 = par[3];
        return getSignalIntegral(mH, width, acc_0, acc_1, acc_2);
    }

    Double_t getTotalIntegral(Double_t mH, Double_t width, Double_t log10_kappa, Double_t acc_0, Double_t acc_1, Double_t acc_2, Double_t a_0, Double_t a_1, Double_t a_2, Double_t a_3, Double_t a_4, Double_t b_0, Double_t b_1, Double_t b_2, Double_t b_3, Double_t b_4)
    {
        Double_t Norm = 6.88734e-08;

        Double_t wfH = width / mH;
        Double_t mh     = 125.09;
        Double_t wfh    = 0.00407 / 125.09;
        Double_t mHc2 = mH*mH*(1+wfH*wfH);
        Double_t mhc2 = mh*mh*(1+wfh*wfh);

        Double_t min_x = mH - 6*width;
        if (min_x<240) min_x=240;
        Double_t max_x = mH + 6*width;
        Double_t n_steps  = 3000;
        Double_t step  = (max_x - min_x) / n_steps;

        Double_t result = 0;
        for (int i_step=0; i_step<n_steps; i_step++)
        {
            Double_t x = min_x + step*(i_step+0.5);
            Double_t x2   = x*x;

            //Propagator part
            // H propagator
            Double_t A_H = 1 + wfH*wfH;
            Double_t B_H = std::pow((x2-mHc2),2) + std::pow((x2*wfH),2);
            Double_t Propagator_H = A_H/B_H;
            // h-H propagator
            Double_t n1 = 1+wfH*wfh;
            Double_t n2 = wfH-wfh;
            Double_t d1 = x2*x2 + x2*x2*wfH*wfh - x2*mhc2 - x2*mHc2 + mhc2*mHc2;
            Double_t d2 = x2*x2*wfH - x2*x2*wfh + x2*mHc2*wfh - x2*mhc2*wfH;
            Double_t Propagator_h_H = 2 * (n1*d1 + n2*d2) / (d1*d1 + d2*d2);
            // Total propagator
            Double_t Propag = Propagator_H + Propagator_h_H * log10_kappa ;

            Double_t G_gg = std::pow(x,3);
            Double_t t_t = std::pow(x/(2*173.34),2);
            Double_t t_b = std::pow(x/(2*4.18),2);
            Double_t A_2 = 0;
            Double_t pi = TMath::Pi();
            if (t_t<=1 && t_b<=1) {
                A_2 = pow( 1/t_t + (t_t-1)/(t_t*t_t)*pow(TMath::ASin(sqrt(t_t)),2) + 1/t_b + (t_b-1)/(t_b*t_b)*pow(TMath::ASin(sqrt(t_b)),2) ,2);
            } else if (t_t<=1 && t_b>1) {
                Double_t log_b = TMath::Log((1+sqrt(1-1/t_b)) / (1-sqrt(1-1/t_b)));
                A_2 = pow( 1/t_t + (t_t-1)/(t_t*t_t)*pow(TMath::ASin(sqrt(t_t)),2) + 1/t_b - 0.25*(t_b-1)/(t_b*t_b)*(log_b*log_b - pi*pi) ,2) + pow( 0.5*(t_b-1)/(t_b*t_b)*pi*log_b ,2);
            } else if (t_t>1 && t_b>1) {
                Double_t log_b = TMath::Log((1+sqrt(1-1/t_b)) / (1-sqrt(1-1/t_b)));
                Double_t log_t = TMath::Log((1+sqrt(1-1/t_t)) / (1-sqrt(1-1/t_t)));
                A_2 = pow( 1/t_b - 0.25*(t_b-1)/(t_b*t_b)*(log_b*log_b - pi*pi) + 1/t_t - 0.25*(t_t-1)/(t_t*t_t)*(log_t*log_t - pi*pi) ,2) + 0.25*pow( (t_b-1)/(t_b*t_b)*pi*log_b + (t_t-1)/(t_t*t_t)*pi*log_t ,2);
            }
            G_gg *= A_2;

            Double_t G_zz = 1;
            G_zz *= std::pow(x,3) * std::pow((1-std::pow(2*91.2/x,2)),1/2) * (1-std::pow(2*91.2/x,2)+0.75*std::pow(2*91.2/x,4)) * (x>2*91.2);

            Double_t Acc =  acc_0 + acc_1*x + acc_2*x*x ;

            Double_t result_tmp =  Norm * L_gg(x) * std::pow(x,1) * G_gg * Propag * G_zz * Acc ;


            // H_B
            Double_t C = 1 / TMath::Power(x,1);
            C = C / ( TMath::Power(x*x - mHc2,2) + TMath::Power(wfH*x*x,2) );
            Double_t C_a = x*x - mHc2 + TMath::Power(wfH*x,2);
            Double_t C_b = wfH * mHc2;
            Double_t A = a_0 + a_1*x + a_2*TMath::Power(x,2) + a_3*TMath::Power(x,3) + a_4*TMath::Power(x,4);
            Double_t B = b_0 + b_1*x + b_2*TMath::Power(x,2) + b_3*TMath::Power(x,3) + b_4*TMath::Power(x,4);

            result_tmp += 10e+07 * L_gg(x) * C * ( C_a * A + C_b * B ) * log10_kappa ;

            if (result_tmp>0) result+= result_tmp;
        }
        result*=step;

        // if (mH==500) std::cout<<"mH = "<<mH<<", width = "<<width<<", log10_kappa = "<<log10_kappa<<", result = "<<result<<std::endl;
        //std::cout<<"[get_stat::getTotalIntegral] Total integral computed: m_H = "<<mH<<" GeV, Width = "<<width<<" GeV, kappa = "<<kappa<<", acc = "<<acc_0<<" "<<acc_1<<" "<<acc_2<<" "<<", Int_a = "<<a_0<<" "<<a_1<<" "<<a_2<<" "<<a_3<<" "<<a_4<<", Int_b = "<<b_0<<" "<<b_1<<" "<<b_2<<" "<<b_3<<" "<<b_4<<" "<<", xs = "<<result<<"\n";
        return result;
    }

    Double_t getTotalIntegralTF2(Double_t *x, Double_t *par)
    {
        // Check the size of x and par
        /*    size_t x_size = sizeof(x)/sizeof(x[0]);
              size_t par_size = sizeof(par)/sizeof(par[0]);
              if( x_size < 3) {
              cerr <<"getTotalIntegralTF3: Not enough dimentions, at least 3, but only"<< x_size <<"given"<<endl;
              return 0;
              }
              if(par_size < 13){
              cerr <<"getTotalIntegralTF3: Not enough parameters, at least 13, but only"<< par_size <<"given"<<endl;
              return 0;
              }
              */    
        Double_t mH = x[0];
        Double_t kappa = x[1];

        Double_t width = mH*par[0];
        Double_t acc_0 = par[1];
        Double_t acc_1 = par[2];
        Double_t acc_2 = par[3];
        Double_t a_0 = par[4];
        Double_t a_1 = par[5];
        Double_t a_2 = par[6];
        Double_t a_3 = par[7];
        Double_t a_4 = par[8];
        Double_t b_0 = par[9];
        Double_t b_1 = par[10];
        Double_t b_2 = par[11];
        Double_t b_3 = par[12];
        Double_t b_4 = par[13];

        return getTotalIntegral(mH, width, kappa, acc_0, acc_1, acc_2,
                a_0, a_1, a_2, a_3, a_4,
                b_0, b_1, b_2, b_3, b_4);
    }

}
