// dllmain.cpp : Defines the entry point for the DLL application.
#include "pch.h"
#include "Eigen/Dense"
#include "nlohmann/json.hpp"
#include "ASRE_Timo_functions.h"
BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}

extern "C" {
    //__declspec(dllexport) int run() {
    __declspec(dllexport) double* run(int nnode, double* meshX, double* meshY, double* meshZ, double* dispV, double* dispL, double* dispT,
        double Eb, double EoverG, double EsNominal, double nis, double dfoot, double bfoot, double ni_foot, double mu_int, double qz_foot
        ) {
        //std::cout << "ASRE run starts" << std::endl;
        VectorXd vecA = VectorXd::Ones(2) * dispV[0];
        VectorXd vecB = VectorXd::Ones(2) * meshX[0];

        VectorXd x_global_foot = VectorXd::Zero(nnode);
        VectorXd y_global_foot = VectorXd::Zero(nnode);
        VectorXd z_global_foot = VectorXd::Zero(nnode);
        VectorXd ucat = VectorXd::Zero(nnode * 6);
        for (int i = 0; i < nnode; i++) {
            x_global_foot(i) = meshX[i];
            y_global_foot(i) = meshY[i];
            z_global_foot(i) = meshZ[i];
            ucat(i * 6) = dispL[i];
            ucat(i * 6 + 1) = dispT[i];
            ucat(i * 6 + 2) = dispV[i];
        }
        double h_el_foot = std::sqrt((x_global_foot(0) - x_global_foot(nnode - 1)) * (x_global_foot(0) - x_global_foot(nnode - 1)) +
            (y_global_foot(0) - y_global_foot(nnode - 1)) * (y_global_foot(0) - y_global_foot(nnode - 1))) / (nnode - 1);
        //std::cout << "x_global_foot: " << std::endl;
        //std::cout << x_global_foot << std::endl;
        //std::cout << "ucat: " << std::endl;
        //std::cout << ucat << std::endl;
        //std::cout << "h_el_foot: " << std::endl;
        //std::cout << h_el_foot << std::endl;
        /*VectorXd z_global_foot = VectorXd::Zero(nnode);
        VectorXd z_global_foot = VectorXd::Zero(nnode);
        VectorXd z_global_foot = VectorXd::Zero(nnode);*/
        /*VectorXd y = VectorXd::Ones(3);*/
        

        int nelement = nnode - 1;
        int dimK = nnode * 6;
        MatrixXd KKfoot = MatrixXd::Zero(dimK, dimK);
        for (int i = 0; i < nelement; i++) {
            double dx = std::sqrt((meshX[i] - meshX[i + 1]) * (meshX[i] - meshX[i + 1]) + (meshY[i] - meshY[i + 1]) * (meshY[i] - meshY[i + 1]));
            MatrixXd KBern3Delt = KBern3D_foot_TIM(Eb, dfoot, bfoot, dx, EoverG, ni_foot);
            KKfoot(seq((i) * 6, (i) * 6 + 11, 1),
                seq((i) * 6, (i) * 6 + 11, 1)) =
                KKfoot(seq((i) * 6, (i) * 6 + 11, 1),
                    seq((i) * 6, (i) * 6 + 11, 1)) + KBern3Delt;
        }
        //std::cout << "KBern3D: " << std::endl;
        double dx = std::sqrt((meshX[0] - meshX[0 + 1]) * (meshX[0] - meshX[0 + 1]) + (meshY[0] - meshY[0 + 1]) * (meshY[0] - meshY[0 + 1]));
        //std::cout << KBern3D_foot_TIM(Eb, dfoot, bfoot, dx, EoverG, ni_foot)<<std::endl;
        
        MatrixXd FLEX_3DOF = calFlexVaziri(z_global_foot, x_global_foot, y_global_foot, EsNominal, nnode, h_el_foot, nis, bfoot);

        //std::cout << "FLEX_3DOF" << std::endl;
        //std::cout << FLEX_3DOF << std::endl;
        double kh_gazetas;
        double kv_gazetas;

        soilsprings_static_foot_gazetas(&kh_gazetas, &kv_gazetas, EsNominal, nis, bfoot, h_el_foot);
        MatrixXd Ks_3DOF = FLEX_3DOF.inverse();
        MatrixXd Ks = MatrixXd::Zero(nnode * 6, nnode * 6);
        MatrixXd FLEX = MatrixXd::Zero(nnode * 6, nnode * 6);
        for (int i = 0; i < nnode; i++) {
            for (int j = 0; j < nnode; j++) {
                Ks(seq(i * 6, 2 + i * 6, 1), seq(j * 6, 2 + j * 6, 1)) =
                    Ks_3DOF(seq(i * 3, 2 + i * 3, 1), seq(j * 3, 2 + j * 3, 1));
                FLEX(seq(i * 6, 2 + i * 6, 1), seq(j * 6, 2 + j * 6, 1)) =
                    FLEX_3DOF(seq(i * 3, 2 + i * 3, 1), seq(j * 3, 2 + j * 3, 1));
            }
            // ADDITIONAL STIFFNESS IN THE DOF 4 TO REMOVE UNCOSTRAINED ROTATION ABOUT THE BEAM AXIS
            Ks(i * 6 + 3, i * 6 + 3) = kv_gazetas * 1000;
            FLEX(i * 6 + 3, i * 6 + 3) = 1.0/(kv_gazetas * 1000.0);
        }
        Ks_3DOF.resize(0, 0);

        //FLEX = Ks.inverse();
        //std::cout << "FLEX after" << std::endl;
        //std::cout << FLEX << std::endl;
        VectorXd P_el = VectorXd::Zero(6 * nnode);
        calc_P_el(&P_el, nnode, qz_foot, h_el_foot);
        //std::cout << "P_el: " << std::endl;
        //std::cout << P_el << std::endl;
        // Because the deadload is unlikely to cause plastic deformation, use elastic solution to estimate initial displacement
        MatrixXd Amatrix = Ks + KKfoot;
        PartialPivLU<MatrixXd> lu = PartialPivLU<MatrixXd>(Amatrix);
        VectorXd u_P_el = lu.solve(P_el);
        Amatrix.resize(0, 0);
        //std::cout << "u_P_el: " << std::endl;
        //std::cout << u_P_el << std::endl;

        //Elastic solution
        //std::cout << "KKfoot: " << std::endl;
        //std::cout << KKfoot << std::endl;
        //std::cout << "Ks: " << std::endl;
        //std::cout << Ks << std::endl;
        VectorXd u_cat_el = lu.solve(P_el + Ks * ucat);
        //std::memcpy(resultArray, u_cat_el.data(), u_cat_el.size() * sizeof(u_cat_el[0]));
        //std::cout << "u_cat_el: " << std::endl;
        //std::cout << u_cat_el << std::endl;
        

        //Prepare for elastioPlaticIteration to solve ucat induced displacement
        VectorXd Kstar = VectorXd::Zero(nnode * 6);
        MatrixXd Lstar = MatrixXd::Zero(nnode * 6, nnode * 6);
        Lstar = FLEX;
        MatrixXd Stiffness = KKfoot;
        /*for (int i = 0; i < nnode * 6; i++) {
            Kstar(i) = 1.0 / FLEX(i, i);
            Lstar(i, i) = 0.0;
            Stiffness(i, i) += Kstar(i);
        }*/
        for (int i = 0; i < nnode; i++) {
            Kstar(i * 6 + 0) = 1.0 / FLEX(i * 6+0, i * 6 + 0);
            Kstar(i * 6 + 1) = 1.0 / FLEX(i * 6 + 1, i * 6 + 1);
            Kstar(i * 6 + 2) = 1.0 / FLEX(i * 6 + 2, i * 6 + 2);
            Kstar(i * 6 + 3) = 1.0 / FLEX(i * 6 + 3, i * 6 + 3);
            Lstar(i * 6 + 0, i * 6 + 0) = 0.0;
            Lstar(i * 6 + 1, i * 6 + 1) = 0.0;
            Lstar(i * 6 + 2, i * 6 + 2) = 0.0;
            Lstar(i * 6 + 3, i * 6 + 3) = 0.0;
            Stiffness(i * 6 + 0, i * 6 + 0) += Kstar(i * 6 + 0);
            Stiffness(i * 6 + 1, i * 6 + 1) += Kstar(i * 6 + 1);
            Stiffness(i * 6 + 2, i * 6 + 2) += Kstar(i * 6 + 2);
            Stiffness(i * 6 + 3, i * 6 + 3) += Kstar(i * 6 + 3);
        }
        /*std::cout << "Kstar: " << std::endl;
        std::cout << Kstar << std::endl;
        std::cout << "Lstar: " << std::endl;
        std::cout << Lstar << std::endl;
        std::cout << "Stiffness: " << std::endl;
        std::cout << Stiffness << std::endl;*/

        VectorXd uinc = u_P_el;
        VectorXd uip = VectorXd::Zero(uinc.size());
        double lim_t_int = 0;
        double lim_c_int = INFINITY;
        /*std::cout << "ucat: " << std::endl;
        std::cout << ucat << std::endl;*/
        elastoPlasticIterationLDLT(Stiffness, KKfoot, Kstar,
            FLEX, Lstar, uinc, ucat, P_el,
            mu_int, lim_t_int, lim_c_int,
            h_el_foot * bfoot,
            uip);
        /*std::cout << "uinc-u_P_el" << std::endl;
        std::cout << uinc-u_P_el << std::endl;*/
        
        VectorXd F_M_deltaT_el_M, F_N_deltaT_el_M, F_S_deltaT_el_M;
        calInternalForces(&F_M_deltaT_el_M, &F_N_deltaT_el_M, &F_S_deltaT_el_M,
            u_P_el, uinc, Eb, EoverG, h_el_foot, dfoot, bfoot, ni_foot, nnode);
        VectorXd epsilon_vector = calculateStrain(&F_S_deltaT_el_M, &F_M_deltaT_el_M,
            &F_N_deltaT_el_M, Eb, EoverG, bfoot, dfoot, ni_foot);





        /*MatrixXd KKextra;
        MatrixXd KKpg = calKKpg(&KKsoil_rid, kv_gazetas, KKfoot, &KKextra, nnode);

        MatrixXd SS = KKfoot + KKextra;
        MatrixXd lamdasts, lamdastd, Kst;
        calLamdastsLamdastdKst(&lamdasts, &lamdastd, &Kst, FLEX, nnode);*/



        //VectorXd result = uinc-u_P_el;
        VectorXd result = epsilon_vector;
        /*std::cout << "result: " << std::endl;
        std::cout << result << std::endl;*/
        //std::memcpy(resultArray, result.data(), result.size() * sizeof(result[0]));
        //resultArray[0] = 0.1415926;


        return result.data();
    }
}

