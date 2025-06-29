#pragma warning(disable : 4996)

#include "udf.h"
#include "mem.h"
#include "sg.h"
#include "math.h"
#include "surf.h"
#include "dpm.h"
#include "stdio.h"

#define liquid_fraction_critical 0.6 
#define t_solidification_rate 0.000312
#define pi 3.1415926
#define casting_velocity 0.0217

double mushy_velocity_narrow(double y)
{
	return (0.02121+(2*(0.64196)/pi)*(0.13262/((4*(pow(y+0.25279, 2)))+pow(0.13262, 2))))+(0.02121+(2*(0.01946)/pi)*(0.13154/((4*(pow(y+0.43874, 2)))+pow(0.13154, 2))))+(0.02121+(2*(-0.64701)/pi)*(0.13547/((4*(pow(y+0.25321, 2)))+pow(0.13547, 2))))-0.042;
}

double mushy_velocity_wide(double x)
{
	return ((2*(0.000880115)/pi)*(0.05199/((4*(pow(x+0.36731, 2)))+pow(0.05199, 2))))+(0.02087+(2*(0.000425913)/pi)*(0.03519/((4*(pow(x+0.3272, 2)))+pow(0.03519, 2))));
}


DEFINE_INIT(shell_setup_m_v, domain)
{
    if (NULLP(user_particle_vars)) 
    {
        Init_User_Particle_Vars(); 
        printf("User particle vars initialized.\n");  
    }

    if (user_particle_vars) 
    {
        strcpy(user_particle_vars[0].name, "shell_zone_flag");
        strcpy(user_particle_vars[0].label, "shell_Zone_Flag");
        strcpy(user_particle_vars[1].name, "shell_zone_time");
        strcpy(user_particle_vars[1].label, "shell_Zone_Time");

        printf("User particle vars set successfully.\n");
    }
    else 
    {
        printf("Error: user_particle_vars is NULL after initialization.\n");
    }
}


DEFINE_DPM_SCALAR_UPDATE(shell_capture_m_v, c, t, initialize, p)
{
    real x_p= P_POS(p)[0];
	real y_p= P_POS(p)[1];
	real z_p= P_POS(p)[2];
    
    double nf_p = 0.65;    
    double wf_ir = 0.1;   
    double wf_or = -0.1;  
    double cap_ind = 5.0;
	
	real v_p_x_direction = P_VEL(p)[0];	
	real v_p_y_direction = P_VEL(p)[1];	
	real v_p_z_direction = P_VEL(p)[2];
	
	real liquid_fraction = C_LIQF(c, t);
	
	real partical_diameter = P_DIAM(p);
	
	real t_time_critical = partical_diameter / t_solidification_rate;
	
	real v_p = sqrt(pow(P_VEL(p)[0], 2) + pow(P_VEL(p)[1], 2) + pow(P_VEL(p)[2], 2));
	
	double m_v_c_narrow = mushy_velocity_narrow(y_p);
	
	double m_v_c_wide = mushy_velocity_wide(y_p);
	
	real v_p_xz = sqrt(pow(P_VEL(p)[0], 2) + pow(P_VEL(p)[2], 2));
	

    if (z_p > 0.63 && y_p > -0.8 && v_p_z_direction > 0)
    {
        if ( liquid_fraction > liquid_fraction_critical && v_p > m_v_c_narrow)
        {
            P_USER_REAL(p, 0) = 0;  
        }
        else
        {
            if ( liquid_fraction < liquid_fraction_critical && v_p < m_v_c_narrow)
            {
                if (P_USER_REAL(p, 0) == 0)
                {
                    P_USER_REAL(p, 0) = 1.0;  
                    P_USER_REAL(p, 1) = P_TIME(p);  
                }

                if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
                {
                    FILE *shellcaptured_mushy_velocity_narrow_1 = fopen("shell_mushy_velocity_narrow_1.log", "a");
                    fprintf(shellcaptured_mushy_velocity_narrow_1, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
                            P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, m_v_c_narrow, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
                    p->stream_index = -1; 
                    MARK_PARTICLE(p, P_FL_REMOVED); 
                    fclose(shellcaptured_mushy_velocity_narrow_1);
                }
            }
        }
	}

	if (z_p > 0.63 && y_p > -0.8 && v_p_z_direction < 0)
	{
		if (liquid_fraction > liquid_fraction_critical && v_p > m_v_c_narrow)
		{
			P_USER_REAL(p, 0) = 0;
		}
		else
		{
			if (liquid_fraction < liquid_fraction_critical && v_p < m_v_c_narrow &&  t_solidification_rate > v_p_xz )
			{
				if (P_USER_REAL(p, 0) == 0)
				{
					P_USER_REAL(p, 0) = 1.0;
					P_USER_REAL(p, 1) = P_TIME(p);
				}

				if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
				{
					FILE *shellcaptured_mushy_velocity_narrow_2 = fopen("shell_mushy_velocity_narrow_2.log", "a");
					fprintf(shellcaptured_mushy_velocity_narrow_2, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
						P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, m_v_c_narrow, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
					p->stream_index = -1;
					MARK_PARTICLE(p, P_FL_REMOVED);
					fclose(shellcaptured_mushy_velocity_narrow_2);
				}
			}
		}

	}
	
	if (z_p > 0.3 && z_p <= 0.63 && y_p > -0.8 && x_p > 0 && v_p_x_direction > 0)
    {
        if ( liquid_fraction > liquid_fraction_critical && v_p > m_v_c_wide)
        {
            P_USER_REAL(p, 0) = 0;  
        }
        else
        {
            if ( liquid_fraction < liquid_fraction_critical && v_p < m_v_c_wide)
            {
                if (P_USER_REAL(p, 0) == 0)
                {
                    P_USER_REAL(p, 0) = 1.0;  
                    P_USER_REAL(p, 1) = P_TIME(p);  
                }

                if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
                {
                    FILE *shellcaptured_mushy_velocity_wide_1 = fopen("shell_mushy_velocity_wide_1.log", "a");
                    fprintf(shellcaptured_mushy_velocity_wide_1, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
                            P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, m_v_c_wide, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
                    p->stream_index = -1; 
                    MARK_PARTICLE(p, P_FL_REMOVED); 
                    fclose(shellcaptured_mushy_velocity_wide_1);
                }
            }
        }
	}
	
	if (z_p > 0.3 && z_p <= 0.63 && y_p > -0.8 && x_p > 0 && v_p_x_direction < 0)
	{
		if (liquid_fraction > liquid_fraction_critical && v_p > m_v_c_wide)
		{
			P_USER_REAL(p, 0) = 0;
		}
		else
		{
			if (liquid_fraction < liquid_fraction_critical && v_p < m_v_c_wide &&  t_solidification_rate > v_p_xz )
			{
				if (P_USER_REAL(p, 0) == 0)
				{
					P_USER_REAL(p, 0) = 1.0;
					P_USER_REAL(p, 1) = P_TIME(p);
				}

				if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
				{
					FILE *shellcaptured_mushy_velocity_wide_2 = fopen("shell_mushy_velocity_wide_2.log", "a");
					fprintf(shellcaptured_mushy_velocity_wide_2, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
						P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, m_v_c_wide, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
					p->stream_index = -1;
					MARK_PARTICLE(p, P_FL_REMOVED);
					fclose(shellcaptured_mushy_velocity_wide_2);
				}
			}
		}

	}

	if (z_p > 0 && z_p <= 0.3 && y_p > -0.8 && x_p > 0 && v_p_x_direction > 0)
    {
        if ( liquid_fraction > liquid_fraction_critical && v_p > casting_velocity)
        {
            P_USER_REAL(p, 0) = 0;  
        }
        else
        {
            if ( liquid_fraction < liquid_fraction_critical && v_p < casting_velocity)
            {
                if (P_USER_REAL(p, 0) == 0)
                {
                    P_USER_REAL(p, 0) = 1.0;  
                    P_USER_REAL(p, 1) = P_TIME(p);  
                }

                if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
                {
                    FILE *shellcaptured_mushy_velocity_wide_3 = fopen("shell_mushy_velocity_wide_3.log", "a");
                    fprintf(shellcaptured_mushy_velocity_wide_3, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
                            P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, casting_velocity, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
                    p->stream_index = -1; 
                    MARK_PARTICLE(p, P_FL_REMOVED); 
                    fclose(shellcaptured_mushy_velocity_wide_3);
                }
            }
        }
	}
	
	if (z_p > 0 && z_p <= 0.3 && y_p > -0.8 && x_p > 0 && v_p_x_direction < 0)
	{
		if (liquid_fraction > liquid_fraction_critical && v_p > casting_velocity)
		{
			P_USER_REAL(p, 0) = 0;
		}
		else
		{
			if (liquid_fraction < liquid_fraction_critical && v_p < casting_velocity &&  t_solidification_rate > v_p_xz )
			{
				if (P_USER_REAL(p, 0) == 0)
				{
					P_USER_REAL(p, 0) = 1.0;
					P_USER_REAL(p, 1) = P_TIME(p);
				}

				if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
				{
					FILE *shellcaptured_mushy_velocity_wide_4 = fopen("shell_mushy_velocity_wide_4.log", "a");
					fprintf(shellcaptured_mushy_velocity_wide_4, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
						P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, casting_velocity, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
					p->stream_index = -1;
					MARK_PARTICLE(p, P_FL_REMOVED);
					fclose(shellcaptured_mushy_velocity_wide_4);
				}
			}
		}

	}

	if (z_p > 0.3 && z_p <= 0.63 && y_p > -0.8 && x_p < 0 && v_p_x_direction < 0)
    {
        if ( liquid_fraction > liquid_fraction_critical && v_p > m_v_c_wide)
        {
            P_USER_REAL(p, 0) = 0;  
        }
        else
        {
            if ( liquid_fraction < liquid_fraction_critical && v_p < m_v_c_wide)
            {
                if (P_USER_REAL(p, 0) == 0)
                {
                    P_USER_REAL(p, 0) = 1.0;  
                    P_USER_REAL(p, 1) = P_TIME(p);  
                }

                if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
                {
                    FILE *shellcaptured_mushy_velocity_wide_5 = fopen("shell_mushy_velocity_wide_5.log", "a");
                    fprintf(shellcaptured_mushy_velocity_wide_5, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
                            P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, m_v_c_wide, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
                    p->stream_index = -1; 
                    MARK_PARTICLE(p, P_FL_REMOVED); 
                    fclose(shellcaptured_mushy_velocity_wide_5);
                }
            }
        }
	}
	
	if (z_p > 0.3 && z_p <= 0.63 && y_p > -0.8 && x_p < 0 && v_p_x_direction > 0)
	{
		if (liquid_fraction > liquid_fraction_critical && v_p > m_v_c_wide)
		{
			P_USER_REAL(p, 0) = 0;
		}
		else
		{
			if (liquid_fraction < liquid_fraction_critical && v_p < m_v_c_wide &&  t_solidification_rate > v_p_xz )
			{
				if (P_USER_REAL(p, 0) == 0)
				{
					P_USER_REAL(p, 0) = 1.0;
					P_USER_REAL(p, 1) = P_TIME(p);
				}

				if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
				{
					FILE *shellcaptured_mushy_velocity_wide_6 = fopen("shell_mushy_velocity_wide_6.log", "a");
					fprintf(shellcaptured_mushy_velocity_wide_6, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
						P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, m_v_c_wide, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
					p->stream_index = -1;
					MARK_PARTICLE(p, P_FL_REMOVED);
					fclose(shellcaptured_mushy_velocity_wide_6);
				}
			}
		}

	}

	if (z_p > 0 && z_p <= 0.3 && y_p > -0.8 && x_p < 0 && v_p_x_direction < 0)
    {
        if ( liquid_fraction > liquid_fraction_critical && v_p > casting_velocity)
        {
            P_USER_REAL(p, 0) = 0;  
        }
        else
        {
            if ( liquid_fraction < liquid_fraction_critical && v_p < casting_velocity)
            {
                if (P_USER_REAL(p, 0) == 0)
                {
                    P_USER_REAL(p, 0) = 1.0;  
                    P_USER_REAL(p, 1) = P_TIME(p);  
                }

                if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
                {
                    FILE *shellcaptured_mushy_velocity_wide_7 = fopen("shell_mushy_velocity_wide_7.log", "a");
                    fprintf(shellcaptured_mushy_velocity_wide_7, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
                            P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, casting_velocity, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
                    p->stream_index = -1; 
                    MARK_PARTICLE(p, P_FL_REMOVED); 
                    fclose(shellcaptured_mushy_velocity_wide_7);
                }
            }
        }
	}
	
	if (z_p > 0 && z_p <= 0.3 && y_p > -0.8 && x_p < 0 && v_p_x_direction > 0)
	{
		if (liquid_fraction > liquid_fraction_critical && v_p > casting_velocity)
		{
			P_USER_REAL(p, 0) = 0;
		}
		else
		{
			if (liquid_fraction < liquid_fraction_critical && v_p < casting_velocity &&  t_solidification_rate > v_p_xz )
			{
				if (P_USER_REAL(p, 0) == 0)
				{
					P_USER_REAL(p, 0) = 1.0;
					P_USER_REAL(p, 1) = P_TIME(p);
				}

				if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
				{
					FILE *shellcaptured_mushy_velocity_wide_8 = fopen("shell_mushy_velocity_wide_8.log", "a");
					fprintf(shellcaptured_mushy_velocity_wide_8, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
						P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, casting_velocity, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
					p->stream_index = -1;
					MARK_PARTICLE(p, P_FL_REMOVED);
					fclose(shellcaptured_mushy_velocity_wide_8);
				}
			}
		}

	}		
	
 if (z_p > 0.63 && y_p < -0.8 && v_p_z_direction > 0)
    {
        if ( liquid_fraction > liquid_fraction_critical && v_p > casting_velocity)
        {
            P_USER_REAL(p, 0) = 0;  
        }
        else
        {
            if ( liquid_fraction < liquid_fraction_critical && v_p < casting_velocity)
            {
                if (P_USER_REAL(p, 0) == 0)
                {
                    P_USER_REAL(p, 0) = 1.0;  
                    P_USER_REAL(p, 1) = P_TIME(p);  
                }

                if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
                {
                    FILE *shellcaptured_mushy_velocity_narrow_3 = fopen("shell_mushy_velocity_narrow_3.log", "a");
                    fprintf(shellcaptured_mushy_velocity_narrow_3, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
                            P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, casting_velocity, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
                    p->stream_index = -1; 
                    MARK_PARTICLE(p, P_FL_REMOVED); 
                    fclose(shellcaptured_mushy_velocity_narrow_3);
                }
            }
        }
	}

	if (z_p > 0.63 && y_p > -0.8 && v_p_z_direction < 0)
	{
		if (liquid_fraction > liquid_fraction_critical && v_p > casting_velocity)
		{
			P_USER_REAL(p, 0) = 0;
		}
		else
		{
			if (liquid_fraction < liquid_fraction_critical && v_p < casting_velocity &&  t_solidification_rate > v_p_xz )
			{
				if (P_USER_REAL(p, 0) == 0)
				{
					P_USER_REAL(p, 0) = 1.0;
					P_USER_REAL(p, 1) = P_TIME(p);
				}

				if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
				{
					FILE *shellcaptured_mushy_velocity_narrow_4 = fopen("shell_mushy_velocity_narrow_4.log", "a");
					fprintf(shellcaptured_mushy_velocity_narrow_4, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
						P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, casting_velocity, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
					p->stream_index = -1;
					MARK_PARTICLE(p, P_FL_REMOVED);
					fclose(shellcaptured_mushy_velocity_narrow_4);
				}
			}
		}

	}

	if (z_p > 0 && y_p < -0.8 && x_p > 0 && v_p_x_direction > 0)
    {
        if ( liquid_fraction > liquid_fraction_critical && v_p > casting_velocity)
        {
            P_USER_REAL(p, 0) = 0;  
        }
        else
        {
            if ( liquid_fraction < liquid_fraction_critical && v_p < casting_velocity)
            {
                if (P_USER_REAL(p, 0) == 0)
                {
                    P_USER_REAL(p, 0) = 1.0;  
                    P_USER_REAL(p, 1) = P_TIME(p);  
                }

                if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
                {
                    FILE *shellcaptured_mushy_velocity_wide_9 = fopen("shell_mushy_velocity_wide_9.log", "a");
                    fprintf(shellcaptured_mushy_velocity_wide_9, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
                            P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, casting_velocity, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
                    p->stream_index = -1; 
                    MARK_PARTICLE(p, P_FL_REMOVED); 
                    fclose(shellcaptured_mushy_velocity_wide_9);
                }
            }
        }
	}
	
	if (z_p > 0 && y_p < -0.8 && x_p > 0 && v_p_x_direction < 0)
	{
		if (liquid_fraction > liquid_fraction_critical && v_p > casting_velocity)
		{
			P_USER_REAL(p, 0) = 0;
		}
		else
		{
			if (liquid_fraction < liquid_fraction_critical && v_p < casting_velocity &&  t_solidification_rate > v_p_xz )
			{
				if (P_USER_REAL(p, 0) == 0)
				{
					P_USER_REAL(p, 0) = 1.0;
					P_USER_REAL(p, 1) = P_TIME(p);
				}

				if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
				{
					FILE *shellcaptured_mushy_velocity_wide_10 = fopen("shell_mushy_velocity_wide_10.log", "a");
					fprintf(shellcaptured_mushy_velocity_wide_10, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
						P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, casting_velocity, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
					p->stream_index = -1;
					MARK_PARTICLE(p, P_FL_REMOVED);
					fclose(shellcaptured_mushy_velocity_wide_10);
				}
			}
		}

	}

	if (z_p > 0 && y_p < -0.8 && x_p < 0 && v_p_x_direction < 0)
    {
        if ( liquid_fraction > liquid_fraction_critical && v_p > casting_velocity)
        {
            P_USER_REAL(p, 0) = 0;  
        }
        else
        {
            if ( liquid_fraction < liquid_fraction_critical && v_p < casting_velocity)
            {
                if (P_USER_REAL(p, 0) == 0)
                {
                    P_USER_REAL(p, 0) = 1.0;  
                    P_USER_REAL(p, 1) = P_TIME(p);  
                }

                if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
                {
                    FILE *shellcaptured_mushy_velocity_wide_11 = fopen("shell_mushy_velocity_wide_11.log", "a");
                    fprintf(shellcaptured_mushy_velocity_wide_11, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
                            P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, casting_velocity, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
                    p->stream_index = -1; 
                    MARK_PARTICLE(p, P_FL_REMOVED); 
                    fclose(shellcaptured_mushy_velocity_wide_11);
                }
            }
        }
	}
	
	if (z_p > 0 && y_p < -0.8 && x_p < 0 && v_p_x_direction > 0)
	{
		if (liquid_fraction > liquid_fraction_critical && v_p > casting_velocity)
		{
			P_USER_REAL(p, 0) = 0;
		}
		else
		{
			if (liquid_fraction < liquid_fraction_critical && v_p < casting_velocity &&  t_solidification_rate > v_p_xz )
			{
				if (P_USER_REAL(p, 0) == 0)
				{
					P_USER_REAL(p, 0) = 1.0;
					P_USER_REAL(p, 1) = P_TIME(p);
				}

				if (P_USER_REAL(p, 0) == 1.0 && P_TIME(p) - P_USER_REAL(p, 1) >= t_time_critical)
				{
					FILE *shellcaptured_mushy_velocity_wide_12 = fopen("shell_mushy_velocity_wide_12.log", "a");
					fprintf(shellcaptured_mushy_velocity_wide_12, "((%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%13.4e\t%9.4f_) injection-0:%d)\n",
						P_POS(p)[0], P_POS(p)[1], P_POS(p)[2], P_VEL(p)[0], P_VEL(p)[1], P_VEL(p)[2], P_DIAM(p), 0.0, casting_velocity, v_p, cap_ind, P_T(p), P_TIME(p), p->part_id);
					p->stream_index = -1;
					MARK_PARTICLE(p, P_FL_REMOVED);
					fclose(shellcaptured_mushy_velocity_wide_12);
				}
			}
		}

	}	
	
}
