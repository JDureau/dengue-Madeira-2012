/**************************************************************************
 *    This file is part of plom.
 *
 *    plom is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published
 *    by the Free Software Foundation, either version 3 of the
 *    License, or (at your option) any later version.
 *
 *    plom is distributed in the hope that it will be useful, but
 *    WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with plom.  If not, see
 *    <http://www.gnu.org/licenses/>.
 *************************************************************************/

#include "simulation.h"

/* automatically generated code: order of the parameters */

#define ORDER_S 0
#define ORDER_I1 1
#define ORDER_I2 2
#define ORDER_R1 3
#define ORDER_R2 4
#define ORDER_S1 5
#define ORDER_S2 6
#define ORDER_I12 7
#define ORDER_I21 8
#define ORDER_beta 9
#define ORDER_psi 10
#define ORDER_i 11
#define ORDER_gamma 12
#define ORDER_alpha 13
#define ORDER_e 14
#define ORDER_d 15
#define ORDER_sto 16
#define ORDER_rep 17
#define ORDER_phi 18


#define ORDER_U 9
#define ORDER_R 10




#define ORDER_N 0
#define ORDER_mu_b 1
#define ORDER_mu_d 2
#define ORDER_prop 3

int step_lyap (double t, const double X[], double f[], void *params)
{
    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_par *p_par = p_calc->p_par;  /* syntaxic shortcut */
    struct s_data *p_data = p_calc->p_data;

    int i;
    int c, ac, cac;    
    double **par = p_par->natural;

    struct s_router **routers = p_data->routers;  /* syntaxic shortcut */

    /* non linear system (automatically generated code)*/
    double _r[N_CAC][6];

    
    double _sf[N_CAC][1];

    

    for(cac=0;cac<N_CAC;cac++) {
	

        
        _sf[cac][0] = sin(2.0*M_PI*(par[ORDER_d][routers[ORDER_d]->map[cac]]+t/ONE_YEAR));

        
        _r[cac][0] = par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)*(X[ORDER_I1*N_CAC+cac]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_I21*N_CAC+cac]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_i][routers[ORDER_i]->map[cac]]);
        _r[cac][1] = gsl_spline_eval(p_calc->spline[ORDER_mu_b][cac],t,p_calc->acc[ORDER_mu_b][cac])*gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _r[cac][2] = par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)*(X[ORDER_I12*N_CAC+cac]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_I2*N_CAC+cac]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_i][routers[ORDER_i]->map[cac]]);
        _r[cac][3] = par[ORDER_alpha][routers[ORDER_alpha]->map[cac]];
        _r[cac][4] = gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac]);
        _r[cac][5] = par[ORDER_gamma][routers[ORDER_gamma]->map[cac]];
    }

    for (c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

	    
	    f[0*N_CAC+cac] =  - (_r[cac][0]*X[ORDER_S*N_CAC+cac]) - (_r[cac][2]*X[ORDER_S*N_CAC+cac]) - (_r[cac][4]*X[ORDER_S*N_CAC+cac]) + (_r[cac][1]);
	    f[1*N_CAC+cac] =  - (_r[cac][5]*X[ORDER_I1*N_CAC+cac]) - (_r[cac][4]*X[ORDER_I1*N_CAC+cac]) + (_r[cac][0]*X[ORDER_S*N_CAC+cac]);
	    f[2*N_CAC+cac] =  - (_r[cac][5]*X[ORDER_I2*N_CAC+cac]) - (_r[cac][4]*X[ORDER_I2*N_CAC+cac]) + (_r[cac][2]*X[ORDER_S*N_CAC+cac]);
	    f[3*N_CAC+cac] =  - (_r[cac][3]*X[ORDER_R1*N_CAC+cac]) - (_r[cac][4]*X[ORDER_R1*N_CAC+cac]) + (_r[cac][5]*X[ORDER_I1*N_CAC+cac]);
	    f[4*N_CAC+cac] =  - (_r[cac][3]*X[ORDER_R2*N_CAC+cac]) - (_r[cac][4]*X[ORDER_R2*N_CAC+cac]) + (_r[cac][5]*X[ORDER_I2*N_CAC+cac]);
	    f[5*N_CAC+cac] =  - (_r[cac][2]*X[ORDER_S1*N_CAC+cac]) - (_r[cac][4]*X[ORDER_S1*N_CAC+cac]) + (_r[cac][3]*X[ORDER_R1*N_CAC+cac]);
	    f[6*N_CAC+cac] =  - (_r[cac][0]*X[ORDER_S2*N_CAC+cac]) - (_r[cac][4]*X[ORDER_S2*N_CAC+cac]) + (_r[cac][3]*X[ORDER_R2*N_CAC+cac]);
	    f[7*N_CAC+cac] =  - (_r[cac][5]*X[ORDER_I12*N_CAC+cac]) - (_r[cac][4]*X[ORDER_I12*N_CAC+cac]) + (_r[cac][2]*X[ORDER_S1*N_CAC+cac]);
	    f[8*N_CAC+cac] =  - (_r[cac][5]*X[ORDER_I21*N_CAC+cac]) - (_r[cac][4]*X[ORDER_I21*N_CAC+cac]) + (_r[cac][0]*X[ORDER_S2*N_CAC+cac]);
        }
    }

    /* linear system: product of jacobian matrix (DIM*DIM) per

       | y[1*DIM+0]       y[1*DIM+1] ...     y[1*DIM+(DIM-1)]   |
       | y[2*DIM+0]       y[2*DIM+1] ...     y[2*DIM+(DIM-1)]   |
       | ...                                                    |
       | y[DIM*DIM+0]     y[DIM*DIM+1] ...  y[DIM*DIM+(DIM-1)]  |

       (automaticaly generated code)
    */

    double _rj[N_CAC][21];
    for(cac=0; cac<N_CAC; cac++){
        
        _rj[cac][0] = X[ORDER_S*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][1] = -X[ORDER_S2*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][2] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_gamma][routers[ORDER_gamma]->map[cac]]+X[ORDER_S2*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][3] = X[ORDER_S1*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][4] = -X[ORDER_S*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][5] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_alpha][routers[ORDER_alpha]->map[cac]];
        _rj[cac][6] = par[ORDER_gamma][routers[ORDER_gamma]->map[cac]];
        _rj[cac][7] = par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)*(X[ORDER_I1*N_CAC+cac]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_I21*N_CAC+cac]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_i][routers[ORDER_i]->map[cac]]);
        _rj[cac][8] = -X[ORDER_S1*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][9] = par[ORDER_alpha][routers[ORDER_alpha]->map[cac]];
        _rj[cac][10] = -X[ORDER_S2*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][11] = X[ORDER_S2*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][12] = -X[ORDER_S1*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][13] = 0;
        _rj[cac][14] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)*(X[ORDER_I1*N_CAC+cac]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_I21*N_CAC+cac]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_i][routers[ORDER_i]->map[cac]]);
        _rj[cac][15] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_gamma][routers[ORDER_gamma]->map[cac]]+X[ORDER_S*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][16] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)*(X[ORDER_I1*N_CAC+cac]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_I21*N_CAC+cac]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_i][routers[ORDER_i]->map[cac]])-par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)*(X[ORDER_I12*N_CAC+cac]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_I2*N_CAC+cac]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_i][routers[ORDER_i]->map[cac]]);
        _rj[cac][17] = -X[ORDER_S*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][18] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_gamma][routers[ORDER_gamma]->map[cac]]+X[ORDER_S1*N_CAC+cac]*par[ORDER_beta][routers[ORDER_beta]->map[cac]]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac]);
        _rj[cac][19] = -gsl_spline_eval(p_calc->spline[ORDER_mu_d][cac],t,p_calc->acc[ORDER_mu_d][cac])-par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)*(X[ORDER_I12*N_CAC+cac]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_I2*N_CAC+cac]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_i][routers[ORDER_i]->map[cac]]);
        _rj[cac][20] = par[ORDER_beta][routers[ORDER_beta]->map[cac]]*(par[ORDER_e][routers[ORDER_e]->map[cac]]*_sf[cac][0]+1)*(X[ORDER_I12*N_CAC+cac]*par[ORDER_psi][routers[ORDER_psi]->map[cac]]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+X[ORDER_I2*N_CAC+cac]/gsl_spline_eval(p_calc->spline[ORDER_N][cac],t,p_calc->acc[ORDER_N][cac])+par[ORDER_i][routers[ORDER_i]->map[cac]]);
    }


    for(c=0; c<N_C; c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;
            for(i=0; i<(N_PAR_SV*N_CAC); i++) {
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (0*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][16])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][17])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][17])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][4])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][4])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (1*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][7])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][15])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][0])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (2*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][20])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][15])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][0])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (3*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][6])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][5])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (4*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][6])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][5])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (5*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][12])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][9])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][19])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][8])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (6*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][10])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][9])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][14])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][1])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (7*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][3])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][20])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][18])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, (8*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = +(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (0*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][11])*X[N_PAR_SV*N_CAC+ (1*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (2*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (3*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (4*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (5*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][7])*X[N_PAR_SV*N_CAC+ (6*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][13])*X[N_PAR_SV*N_CAC+ (7*N_CAC+ cac)*N_PAR_SV*N_CAC +i]+(_rj[cac][2])*X[N_PAR_SV*N_CAC+ (8*N_CAC+ cac)*N_PAR_SV*N_CAC +i];
                
            }
        }
    }

    return GSL_SUCCESS;
}
