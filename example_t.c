
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eph_manager.h"
#include "novas.h"
void juliandaytoutc(double *ju, short int *year, short int *month, short int *day, short int *hour, short int *min, double *second);

int main (void)
{
    
    // UTC date and time
    short int year = 2019; //UTC
    short int month = 7;
    short int day = 2;
	double hour = 15.0;  //hour
	short int hour1;
    const short int leap_secs = 37; //sec 1921=0.0
	short int min = 0;
	double sec = 0.0;
    //EOP
    const double ut1_utc = 0.06723; /* UT1-UTC (s)  */
    const double x_pole = 0.1538; /* polar motion (arcsec) */
    const double y_pole = 0.4389;

    //accuracy
   const short int accuracy = 0; //0 ... full accuracy, 1 ... reduced accuracy equ2hor
   short int error = 0;
   short int de_num = 0;
    short int coord_sys=1;  //0 ... GCRS or "local GCRS", 1 ...True equator and equinox of date, = 3 ... astrometric coordinates, i.e., without light deflection or aberration.
    
    //大金山岛(最高点)lat= 121.4210378600   lon=  30.6915202617  height= 103.7000000000
    //佘山岛（崇明南）lat= 122.2408852780   lon=  31.4216961532  height=   2.1900000000
    //南汇嘴(自贸区),121.97214203500,30.88479859870，2.19
    //C919起飞地(浦东机场第四跑道),121.82640451500,31.13313227050，2.19
    //上海中心,121.50108770200,31.23600120130，632.0，2.19
    //海关大钟,121.48549131700,31.23875973740，79.2
    //一大会址,121.47072041000,31.22206280900，2.19
    //四叶草,121.30235221000,31.19140111610，2.19
    //江浙沪交界,120.86691019000,31.12358387180，2.19
    //上海光源,121.58411631200,31.19062470030，2.19
    //佘山天文台,121.186071688954,31.0959623108213，99.26
    //鸡骨礁（上海全境最东端的岛礁）东经122.3825277777778 122°22′57.1，北纬31.174083333333332 31°10′26.7
    
    double ddzz=0.00084; //deg 0.00084
    double ddsec=0.001; //sec 0.001
 
    // observer on the surface of the Earth
   const double latitude = 31.174083333333332;  /* latitude (degree) */
   const double longitude = 122.3825277777778; /* east longitude (degree) */
   const double height = 12.2; /* height a.s.l (m) */
   const double temperature = 25.0; /* temperature T (K=273.15+T) */
   const double pressure = 1010.0; /* pressure (hPa) */
   const short int ref_option = 2; /*0 ... no refraction, 1 ... include refraction, using 'standard' atmospheric
    conditions, 2 ... include refraction, using atmospheric parameters input in the 'location' structure. */

    //observer on the near-Earth spacecraft1. Both input vectors are with respect to true equator and equinox of date.
   double rsun=696000.0, ddrsun=0.0;
   double rearth=6371.0084; //km
   double he=height/1000.0;
   double dtheta=0.0;
   //const double x_pole = -0.002;
   //const double y_pole = +0.529;
//   double jd_beg, jd_end, jd_utc, jd_tt, jd_ut1, x, secdif, jd_tdb, delta_t, zd, az, rar, decr;
   double jd_beg, jd_end, jd_utc, jd_tt, jd_ut1, x, secdif, jd_tdb, delta_t;
   //structures
   observer obs_geoc, obs_loc ;
   cat_entry cat_star, dummy_star;
   object star, moon, mars, neptune, sun,earth;
   sky_pos t_place;
/* 计算日月食迭代要用到的变量 */
   static double peb[3], veb[3], peb0[3], veb0[3], pmb[3], vmb[3], jd[2], jd0[2];//jd代表当前时间，peb和veb、pmb和vmb代表jd时间的地球或者月亮相对于太阳的位置和速度矢量，
   //jd0代表光离开太阳的时间（约jd前8分钟），peb0和veb0代表jd0时间的地球或者月亮相对于太阳的位置和速度矢量，
   double rS,rE, rM, TA;
   double theta, theta_P, thetaE, thetaM, OElength, OMlength, thetaE_P, thetaM_P, OElength_P, OMlength_P, SElength, SMlength;
   double R0[3], R0P[3], OE[3], OEP[3], OM[3], MO[3], ME[3], OMP[3], OE_S[3], OE_ST[3], OT[3], OTT[3], ET[3], ETT[3], SE[3]; //对应ppt图上的点
   double vpd, vpd_P, vpd_ME;   //向量积
   //double theta1, s, OU[3], OR[3];
//   double moon_sun_dis;
   TA = AU / C;
   //距离以AU为单位
   rE = (ERAD + 65000.0 )/ AU;//考虑地球大气层的影响
   rE = (ERAD) / AU;
   rS = 696000000.0 / AU;
   rM = 1738000.0 / AU;
/* Make a structure of the observer */
   make_observer_at_geocenter ( &obs_geoc );
   make_observer_on_surface (latitude, longitude, height, temperature, pressure, &obs_loc);
    
/* Make structures of type 'object' for the Star, Sun, Moon, Mars and earth. */
    make_cat_entry ("647080","SC10MA",647080,0.001823085,25.88645705, 20.23, -7.14, 4.34, -31.0, &cat_star);
    make_cat_entry ("DUMMY","xxx",0,0.0,0.0,0.0,0.0,0.0,0.0,&dummy_star);
    
    if ((error = make_object (2,0,"star",&cat_star, &star)) != 0)
    {
        printf ("Error %d from make_object (cat_star)\n", error);
        return (error);
    }
    if ((error = make_object (0,11,"Moon",&dummy_star, &moon)) != 0) // Mercury = 1,...,Pluto = 9, Sun = 10, Moon = 11
   {
      printf ("Error %d from make_object (Moon)\n", error);
      return (error);
   }
	if ((error = make_object(0, 3, "Earth", &dummy_star, &earth)) != 0) // Mercury = 1,...,Pluto = 9, Sun = 10, Moon = 11
	{
		printf("Error %d from make_object (Earth)\n", error);
		return (error);
	}
   if ((error = make_object (0,4,"Mars",&dummy_star, &mars)) != 0)
   {
      printf ("Error %d from make_object (Mars)\n", error);
      return (error);
   }
    if ((error = make_object (0,8,"Neptune",&dummy_star, &neptune)) != 0)
    {
        printf ("Error %d from make_object (Neptune)\n", error);
        return (error);
    }
    if ((error = make_object (0,10,"Sun",&dummy_star, &sun)) != 0)
    {
        printf ("Error %d from make_object (Sun)\n", error);
        return (error);
    }
/* Open the JPL binary ephemeris file, here named "JPLEPH". */
   if ((error = ephem_open ("JPLEPH", &jd_beg,&jd_end,&de_num)) != 0)
   {
      if (error == 1)
         printf ("JPL ephemeris file not found.\n");
       else
         printf ("Error reading JPL ephemeris file header.\n");
      return (error);
   }
    else
   {
      //printf ("JPL ephemeris DE%d open. Start JD = %10.2f  End JD = %10.2f\n",
       //  de_num, jd_beg, jd_end);
     // printf ("\n");
   }

    hour = hour + ddsec/3600.0 ;
    
    /* Establish time arguments. */
    jd_utc = julian_date (year,month,day,hour);
	//printf("jd_utc:%f\n", jd_utc);
	double s1 = 10.0, s2 = 0.0;
//	double elon, elat;
	jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
	jd_ut1 = jd_utc + ut1_utc / 86400.0;
	delta_t = 32.184 + leap_secs - ut1_utc;
	secdif = 0.0;
	tdb2tt(jd_tt, &x, &secdif);
	jd_tdb = jd_tt + secdif / 86400.0;
	jd[0] = jd_tdb;
	jd[1] = 0.0;
	jd0[0] = jd[0] - TA / 86400.0;
	jd0[1] = 0.0;
//24节气
	/*do{
		jd_utc += 0.00001;
		jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
		delta_t = 32.184 + leap_secs - ut1_utc;
		if ((error = place(jd_tt, &sun, &obs_geoc, delta_t, coord_sys, accuracy, &t_place)) != 0)
		{
			printf("Error %d from place.", error);
			return (error);
		}
		equ2ecl(jd_tt, coord_sys, accuracy, t_place.ra, t_place.dec, &elon, &elat);
		//printf("黄经:%5.2f，黄纬:%5.2f\n", elon,elat);
		//printf("fabs(elon-195):%5.2f\n", fabs(elon - 195));
	} while (fabs(elon-225.0) > 0.00001);
	short int hour2;
	jd_utc += 8.0 / 24.0;
	juliandaytoutc(&jd_utc, &year, &month, &day, &hour2, &min, &sec);
	printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n                            2018年立冬时间：%d-%d-%d %d:%02d:%08.6f\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n", year, month, day, hour2, min, sec);
	return;*/

//日食
	//日偏食
	int flag1 = 0, flag2 = 0, flag3 = 0,flag4=0;
	double t1,dt_theta=TWOPI,lon,lat,ratio;
	while (1)
	{
		//printf("s1：%f\n",s1);
		//printf("s2：%f\n", s2);
		jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
		jd_ut1 = jd_utc + ut1_utc / 86400.0;
		delta_t = 32.184 + leap_secs - ut1_utc;
		secdif = 0.0;
		tdb2tt(jd_tt, &x, &secdif);
		jd_tdb = jd_tt + secdif / 86400.0;
		jd[0] = jd_tdb;
		jd[1] = 0.0;
		if ((error = ephemeris(jd, &moon, 1, accuracy, peb, veb)) != 0)
		{
			printf("Error %d from ephemeris_moon", error);
			return (error);
		}
		jd0[0] = jd[0] - TA / 86400.0;
		jd0[1] = 0.0;
		if ((error = ephemeris(jd0, &moon, 1, accuracy, peb0, veb0)) != 0)
		{
			printf("Error %d from ephemeris_moon", error);
			return (error);
		}
		if ((error = ephemeris(jd, &earth, 1, accuracy, pmb, vmb)) != 0)
		{
			printf("Error %d from ephemeris_earth", error);
			return (error);
		}
		for (int i = 0; i < 3; i++)
		{
			//偏食
			R0P[i] = 1 / (1 + rS / rM) * peb0[i] - peb[i]; //最重要一步，RO代表月球半影交点（在太阳和月亮之间）指向太阳质心的矢量
			OMP[i] = R0P[i] + peb[i];//月球半影交点O指向月球质心的矢量，即是1 / (1+rS/rM) * peb0[i]
			//全食或者环食
			R0[i] = 1 / (rM / rS - 1) * peb0[i] + peb0[i] - peb[i]; ////最重要一步，RO代表月球本影交点（在太阳和月亮连线之外）指向太阳质心的矢量
			OM[i] = R0[i] + peb[i];//月球全影交点O指向月球质心的矢量
			MO[i] = -OM[i];
			//月球半影影锥锥点和全影影锥锥点指向地心的矢量
			OEP[i] = R0P[i] + pmb[i]; //偏食月球半影影锥锥点到地心矢量
			OE[i] = R0[i] + pmb[i];//全食月球全影影锥锥点到地心矢量
			ME[i] = MO[i] + OE[i];
			//地球全影影锥锥点指向地心的矢量
			//OE_S[i] = 1 / (rE / rS - 1) * pmb[i];//和月食时候一样
			SE[i] = OEP[i] + peb0[i] - OMP[i];//8分钟前太阳质心到此刻地心的矢量
			OE_S[i] = -1 / (rS / rE - 1) * SE[i];//偏食地球全影影锥锥点到地心矢量
			OE_ST[i] = 1 / (rS / rE + 1) * SE[i];//全食地球半影影锥锥点到地心矢量
			OT[i] = OE_S[i] - OEP[i];//同时切过太阳、月亮、地球边缘的矢量，月球半影影锥和地球全影影锥重合
			OTT[i] = OE_ST[i] - OE[i];//同时切过太阳、月亮、地球边缘的矢量，月球全影影锥和地球半影影锥重合
		}
		//全食
		OElength = sqrt(OE[1] * OE[1] + OE[2] * OE[2] + OE[0] * OE[0]);
		OMlength = sqrt(OM[1] * OM[1] + OM[2] * OM[2] + OM[0] * OM[0]);
		thetaE = asin(rE / OElength);
		thetaM = asin(rM / OMlength);
		//偏食
		OElength_P = sqrt(OEP[1] * OEP[1] + OEP[2] * OEP[2] + OEP[0] * OEP[0]);
		OMlength_P = sqrt(OMP[1] * OMP[1] + OMP[2] * OMP[2] + OMP[0] * OMP[0]);
		thetaE_P= asin(rE / OElength_P);
		thetaM_P = asin(rM / OMlength_P);
		//太阳月球联系和太阳地球连线
		SElength = sqrt(SE[1] * SE[1] + SE[2] * SE[2] + SE[0] * SE[0]);
		SMlength = sqrt(peb0[1] * peb0[1] + peb0[2] * peb0[2] + peb0[0] * peb0[0]);
		//SElength = sqrt(ME[1] * ME[1] + ME[2] * ME[2] + ME[0] * ME[0]);
		//SMlength = sqrt(MO[1] * MO[1] + MO[2] * MO[2] + MO[0] * MO[0]);
		vpd = 0.0;
		vpd_P = 0.0;
		vpd_ME = 0.0;
		for (int i = 0; i < 3; i++)
		{
			vpd += OM[i] * OE[i];
			vpd_P += OMP[i] * OEP[i];
			vpd_ME += peb0[i] * SE[i];
			//vpd_ME += ME[i] * MO[i];
		}
		//偏食
		theta_P = acos(vpd_P / (OElength_P*OMlength_P));
		if (OElength_P > OMlength_P &&fabs((theta_P - (thetaE_P + thetaM_P))) <= 1e-6&&fabs(t1 - jd_utc)>60.0 / 86400.0)//开始结束时间差大于1分钟
		{
			//printf("半影影锥长度：%f\n", rM / OMlength_P*(OMlength_P+380000000/AU)*AU);
			//求月球半影影锥侧线刚切到地球表面时的点的地理坐标
			ratio = (OE_S[0] * OT[0] + OE_S[1] * OT[1] + OE_S[2] * OT[2]) / (OT[0] * OT[0] + OT[1] * OT[1] + OT[2] * OT[2]);
			for (int i = 0; i < 3; i++)
			{
				ET[i] = ratio*OT[i] - OE_S[i];
			}
			cel2ter(jd_ut1, 0.0, delta_t, 1, accuracy, 1, x_pole, y_pole, ET, ETT);
			vector2radec(ETT, &lon, &lat);
			//月球位于正天顶的地球坐标
			if ((error = place(jd_tt, &moon, &obs_geoc, delta_t, coord_sys, accuracy, &t_place)) != 0)//topocentric place: location->where: on surface of earth , coord_sys : true equator and equinox of date
			{
				printf("Error %d from place.", error);
				return (error);
			}
			if (lon > 12.0)lon=lon-24.0;
			t1 = jd_utc; 
			jd_utc += 8.0 / 24.0;
			juliandaytoutc(&jd_utc, &year, &month, &day, &hour1, &min, &sec);
			if (flag1 == 1)
			{ 
				printf("%d年%d月%d日日偏食结束的北京时间： %02d:%02d:%02.0f,半影锥表面切到地球表面的坐标:经度：%f,纬度：%f\n", year, month, day, hour1, min, sec,lon*15.0,lat);
				break;
			}
			else 
			{
				printf("%d年%d月%d日日偏食开始的北京时间： %02d:%02d:%02.0f,半影锥表面切到地球表面的坐标:经度：%f,纬度：%f\n", year, month, day, hour1, min, sec, lon*15.0, lat);
			};
			jd_utc -= 8.0 / 24.0;
			flag1++;
		}
		//全食
		theta = acos(vpd / (OElength*OMlength));
		if (OElength < OMlength &&fabs((theta - (thetaE + thetaM))) <= 1e-5 && fabs(t1 - jd_utc)>60.0 / 86400.0)
		{
			//求月球全影影锥侧线刚切到地球表面时的点的地理坐标
			ratio = (OE_ST[0] * OTT[0] + OE_ST[1] * OTT[1] + OE_ST[2] * OTT[2]) / (OTT[0] * OTT[0] + OTT[1] * OTT[1] + OTT[2] * OTT[2]);
			for (int i = 0; i < 3; i++)
			{
				ET[i] = ratio*OTT[i] - OE_ST[i];
			}
			cel2ter(jd_ut1, 0.0, delta_t, 1, accuracy, 1, x_pole, y_pole, ET, ETT);
			vector2radec(ETT, &lon, &lat);
			if (lon > 12.0)lon = lon - 24.0;
			t1 = jd_utc; 
			jd_utc += 8.0 / 24.0;
			juliandaytoutc(&jd_utc, &year, &month, &day, &hour1, &min, &sec);
			if (flag2 == 1)
			{
				printf("%d年%d月%d日日全食结束的北京时间： %02d:%02d:%02.0f,全影锥表面切到地球表面的坐标:经度：%f,纬度：%f\n", year, month, day, hour1, min, sec, lon*15.0, lat);
			}
			else{ printf("%d年%d月%d日日全食结束的北京时间： %02d:%02d:%02.0f,全影锥表面切到地球表面的坐标:经度：%f,纬度：%f\n", year, month, day, hour1, min, sec, lon*15.0, lat); };
			jd_utc -= 8.0 / 24.0;
			flag2++;
		}
		if (OElength < OMlength &&fabs(TWOPI / 2.0 - (theta + thetaE + thetaM)) <= 1e-5 && fabs(t1 - jd_utc)>60.0 / 86400.0){
			//dt_theta = 0.0;//环食食甚的时候OE和OR夹角最大，设置初始值足够小
			t1 = jd_utc;
			jd_utc += 8.0 / 24.0;
			juliandaytoutc(&jd_utc, &year, &month, &day, &hour1, &min, &sec);
			if (flag3 == 1)
			{
				printf("%d年%d月%d日日环食结束的北京时间： %02d:%02d:%02.0f\n", year, month, day, hour1, min, sec);
			}
			else{ printf("%d年%d月%d日日环食开始的北京时间： %02d:%02d:%02.0f\n", year, month, day, hour1, min, sec); };
			jd_utc -= 8.0 / 24.0;
			flag3++;
		}
		//printf("theta：%20.19f\n", theta);
		//printf("dt_theta：%20.19f\n", dt_theta);
		if (dt_theta<theta_P&&flag4==0){
			ratio = (vpd_ME - sqrt(vpd_ME*vpd_ME - SMlength *SMlength*(SElength*SElength- rE*rE))) / (SMlength *SMlength);
			printf("ratio:%6.4f\n", ratio);
			for (int i = 0; i < 3; i++)
			{
				ET[i] = ratio*peb0[i] - SE[i];
			}
			cel2ter(jd_ut1, 0.0, delta_t, 1, accuracy, 1, x_pole, y_pole, ET, ETT);
			vector2radec(ETT, &lon, &lat);
			if (lon > 12.0)lon = lon - 24.0;
			if ((error = place(jd_tt, &moon, &obs_geoc, delta_t, coord_sys, accuracy, &t_place)) != 0)//topocentric place: location->where: on surface of earth , coord_sys : true equator and equinox of date
			{
				printf("Error %d from place.", error);
				return (error);
			}
			//sidereal_time(jd_ut1, 0.0, delta_t, 1, 1, 1, &lon);
			//lon = lon - t_place.ra;
			//if (t_place.ra > 12.0)lon = t_place.ra - 24.0;
			jd_utc += 8.0 / 24.0;
			juliandaytoutc(&jd_utc, &year, &month, &day, &hour1, &min, &sec);
			printf("%d年%d月%d日食甚的北京时间： %02d:%02d:%02.0f,锥线在地球投影坐标:经度：%f,纬度：%f\n", year, month, day, hour1, min, sec,lon*15.0,lat);
			jd_utc -= 8.0 / 24.0;
			flag4++;
		}
		dt_theta = theta_P;
		jd_utc += 0.000001;
		//printf("vpd：%10.9f\n", vpd);
		//printf("OElength：%10.9f\n", OElength);
		//printf("OMlength：%10.9f\n", OMlength);
		//printf("OElength*OMlength：%10.9f\n", (OElength*OMlength));
		//printf("vpd / OElength*OMlength：%10.9f\n", vpd / (OElength*OMlength));
		//theta = acos(vpd / (OElength*OMlength));
		//printf("theta：%10.9f\n", theta);
		//s2 = theta ;
		//if (s2 > s1)break;
		//short int hour2;
		//jd_utc += 8.0 / 24.0;
		//juliandaytoutc(&jd_utc, &year, &month, &day, &hour2, &min, &sec);
		//printf("时间：%f\n", jd_tdb);
		//printf("时间：%d-%d-%d-%d:%02d:%02.2f\n", year, month, day, hour2, min, sec);
		//printf("OElength：%f\n", OElength);
		//printf("OMlength：%f\n", OMlength);
		//printf("vpd：%f\n", vpd);
		//jd_utc += 0.0001;
		//printf("jd[0]：%f\n", jd[0]);
		//printf("theta：%20.19f\n", theta);
		//printf("thetaE：%20.19f\n", thetaE);
		//printf("thetaM：%20.19f\n",thetaM);
		//printf("deltseta：%20.19f\n", fabs((theta - (thetaE + thetaM))));
		//jd_utc -= 8.0 / 24.0;
		//s1 = s2;
	//} while (OElength < OMlength || fabs((theta - (thetaE + thetaM))) >= 1e-6);
	}// while (OElength > OMlength || fabs((theta - (thetaE + thetaM))) >= 1e-5);
	printf("##########################################end###############################################\n");
	//jd_utc += 8.0 / 24.0;
	//short int hour1;
	//juliandaytoutc(&jd_utc, &year, &month, &day, &hour1, &min, &sec);
	//printf("日偏食时间：%d-%d-%d %d:%02d:%02.0f\n", year, month, day, hour1, min, sec);
	//jd_utc -= 8.0 / 24.0;
	//日全食或者日环食

	/*double r1, r2;//日心到锥点和地球的距离
	double dtheta1;
	do{
		//printf("s1：%f\n",s1);
		//printf("s2：%f\n", s2);
		jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
		jd_ut1 = jd_utc + ut1_utc / 86400.0;
		delta_t = 32.184 + leap_secs - ut1_utc;
		secdif = 0.0;
		tdb2tt(jd_tt, &x, &secdif);
		jd_tdb = jd_tt + secdif / 86400.0;
		jd[0] = jd_tdb;
		jd[1] = 0.0;
		if ((error = ephemeris(jd, &moon, 1, accuracy, peb, veb)) != 0)
		{
			printf("Error %d from ephemeris_moon", error);
			return (error);
		}
		jd0[0] = jd[0] - TA / 86400.0;
		jd0[1] = 0.0;
		if ((error = ephemeris(jd0, &moon, 1, accuracy, peb0, veb0)) != 0)
		{
			printf("Error %d from ephemeris_moon", error);
			return (error);
		}
		r1 = 0.0;
		r2 = 0.0;
		for (int i = 0; i < 3; i++)
		{
			R0[i] = 1 / (rM / rS - 1) * peb0[i] + peb0[i] - peb[i]; ////最重要一步，RO代表月球本影交点（在太阳和月亮连线之外）指向太阳质心的矢量
			OM[i] = R0[i] + peb[i];
			r1 += R0[i] * R0[i];
		}
		if ((error = ephemeris(jd, &earth, 1, accuracy, pmb, vmb)) != 0)
		{
			printf("Error %d from ephemeris_earth", error);
			return (error);
		}
		for (int i = 0; i < 3; i++)
		{
			OE[i] = R0[i] + pmb[i];
			r2 += pmb[i] * pmb[i];
			//printf("pmb[i]:%f\n", pmb[i]);
		}
		OElength = sqrt(OE[1] * OE[1] + OE[2] * OE[2] + OE[0] * OE[0]);
		OMlength = sqrt(OM[1] * OM[1] + OM[2] * OM[2] + OM[0] * OM[0]);
		thetaE = asin(rE / OElength);
		thetaM = asin(rM / OMlength);
		vpd = 0.0;
		for (int i = 0; i < 3; i++)
		{
			vpd += OM[i] * OE[i];
		}
		theta = acos(vpd / (OElength*OMlength));
		s2 = theta;
		if (r1 < r2)//太阳到锥点的距离小于到地心的距离，日全食开始时OM和OE之间夹角加上地球和月亮视半径角距为180度
		{
			dtheta1 = fabs(TWOPI / 2 - theta - thetaE - thetaM);
		}
		else//太阳到锥点的距离小于到地心的距离，日全食开始时OM和OE之间夹角等于地球和月亮视半径角距之和
		{
			//printf("r1:%f r2:%f\n", r1, r2);
			dtheta1 = fabs((theta - (thetaE + thetaM)));
		}
		jd_utc += 0.000001;
		short int hour2;
		jd_utc += 8.0 / 24.0;
		juliandaytoutc(&jd_utc, &year, &month, &day, &hour2, &min, &sec);
		jd_utc -= 8.0 / 24.0;
		s1 = s2;
	} while (OElength > OMlength || dtheta1 >= 0.00001);
	//} while (OElength > OMlength || dtheta1 >= 0.00001);
	//jd_utc += 8.0 / 24.0;
	juliandaytoutc(&jd_utc, &year, &month, &day, &hour1, &min, &sec);
	//printf("日全食时间：%d-%d-%d %d:%02d:%02.0f\n", year, month, day, hour1, min, sec);
	*/
	//计算月食过程
	//theta = TWOPI;//先使theta1的初值足够大
	/*do{
			if ((error = ephemeris(jd0, &earth, 1, accuracy, peb0, veb0)) != 0)
			{
				printf("Error %d from ephemeris_earth", error);
				return (error);
			}
			if ((error = ephemeris(jd, &earth, 1, accuracy, peb, veb)) != 0)
			{
				printf("Error %d from ephemeris_earth", error);
				return (error);
			}
			for (int i = 0; i < 3; i++)
			{
				R0[i] = 1 / (rE / rS - 1) * peb0[i] + peb0[i] - peb[i]; //R1 - R2是偏差, OE[i] =1 / (rE / rS - 1) * peb0[i],相似三角形
				OE[i] = R0[i] + peb[i];
			}
			if ((error = ephemeris(jd, &moon, 1, accuracy, pmb, vmb)) != 0)
			{
				printf("Error %d from ephemeris_moon", error);
				return (error);
			}
			for (int i = 0; i < 3; i++)
			{
				OM[i] = R0[i] + pmb[i];
			}
			vpd = 0.0;
			for (int i = 0; i < 3; i++)
			{
				vpd += OM[i] * OE[i];
			}
			OElength = sqrt(OE[1] * OE[1] + OE[2] * OE[2] + OE[0] * OE[0]);
			OMlength = sqrt(OM[1] * OM[1] + OM[2] * OM[2] + OM[0] * OM[0]);
			thetaE = asin(rE / OElength);
			thetaM = asin(rM / OMlength);
			//printf("OElength：%f\n", OElength);
			//printf("OMlength：%f\n", OMlength);
			//printf("vpd ：%30.29f\n", vpd);
			theta = acos(vpd / (OElength*OMlength));
			jd_utc += 0.000001;
			//jd_utc += 0.0001;
			jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
			jd_ut1 = jd_utc + ut1_utc / 86400.0;
			delta_t = 32.184 + leap_secs - ut1_utc;
			secdif = 0.0;
			tdb2tt(jd_tt, &x, &secdif);
			jd_tdb = jd_tt + secdif / 86400.0;
			jd[0] = jd_tdb;
			jd[1] = 0.0;
			jd0[0] = jd[0] - TA / 86400.0;
			jd0[1] = 0.0;
			//printf("jd[0]：%f jd0[0]：%f\n", jd[0], jd0[0]);
			//printf("thetaE：%f\n", thetaE);
			//printf("thetaM：%f\n",thetaM);
			//printf("deltseta：%f\n", fabs((theta - (thetaE + thetaM))));
	} while(OElength < OMlength || fabs((theta - (thetaE + thetaM))) >= 0.00001);//月全食迭代条件*/
	//计算满月过程
	/*theta = TWOPI;//先使theta1的初值足够大
	do{
		theta1 = theta;
		if ((error = ephemeris(jd0, &earth, 1, accuracy, peb0, veb0)) != 0)
		{
			printf("Error %d from ephemeris_earth", error);
			return (error);
		}
		if ((error = ephemeris(jd, &earth, 1, accuracy, peb, veb)) != 0)
		{
			printf("Error %d from ephemeris_earth", error);
			return (error);
		}
		for (int i = 0; i < 3; i++)
		{
			R0[i] = 1 / (rE / rS - 1) * peb0[i] + peb0[i] - peb[i]; //R1 - R2是偏差, OE[i] =1 / (rE / rS - 1) * peb0[i],相似三角形
			OE[i] = R0[i] + peb[i];
		}
		//OU=OExR0
		OU[0] = OE[1] * R0[2] - OE[2] * R0[1];
		OU[1] = OE[2] * R0[0] - OE[0] * R0[2];
		OU[2] = OE[0] * R0[1] - OE[1] * R0[0];
		//OR=OExOU
		OR[0] = OE[1] * OU[2] - OE[2] * OU[1];
		OR[1] = OE[2] * OU[0] - OE[0] * OU[2];
		OR[2] = OE[0] * OU[1] - OE[1] * OU[0];
		s = sqrt(OR[1] * OR[1] + OR[2] * OR[2] + OR[0] * OR[0]);
		OR[0] = OR[0] / s;
		OR[1] = OR[1] / s;
		OR[2] = OR[2] / s;
		if ((error = ephemeris(jd, &moon, 1, accuracy, pmb, vmb)) != 0)
		{
			printf("Error %d from ephemeris_moon", error);
			return (error);
		}
		for (int i = 0; i < 3; i++)
		{
			OM[i] = R0[i] + pmb[i];
		}
		s = sqrt(OM[1] * OM[1] + OM[2] * OM[2] + OM[0] * OM[0]);
		OM[0] = OM[0] / s;
		OM[1] = OM[1] / s;
		OM[2] = OM[2] / s;
		vpd = 0.0;
		for (int i = 0; i < 3; i++)
		{
			//printf("peb0%d ：%20.19f\n", i, peb0[i]);
			//printf("peb%d ：%20.19f\n", i, peb[i]);
			//printf("OU%d ：%20.19f\n", i, OU[i]);
			//printf("OR%d ：%20.19f\n", i, OR[i]);
			//printf("vpd%d ：%20.19f\n",i, OE[i] * OR[i]);
			vpd += OM[i] * OR[i];
		}
		OElength = sqrt(OE[1] * OE[1] + OE[2] * OE[2] + OE[0] * OE[0]);
		OMlength = sqrt(OM[1] * OM[1] + OM[2] * OM[2] + OM[0] * OM[0]);
		thetaE = asin(rE / OElength);
		thetaM = asin(rM / OMlength);
		//vpd = 0.0;
		//for (int i = 0; i < 3; i++)
		//{
		//vpd+= OM[i]*OE[i];
		//}
		//OElength = sqrt(peb0[1] * peb0[1] + peb0[2] * peb0[2] + peb0[0] * peb0[0]);
		//OMlength = sqrt(peb[1] * peb[1] + peb[2] * peb[2] + peb[0] * peb[0]);
		//vpd = 0.0;
		//for (int i = 0; i < 3; i++)
		//{
		//vpd += peb0[i] *peb[i];
		//}
		//printf("OElength：%f\n", OElength);
		//printf("OMlength：%f\n", OMlength);
		//printf("vpd ：%30.29f\n", vpd);
		theta = acos(vpd / (OElength*OMlength));
		jd_utc += 0.000001;
		//jd_utc += 0.0001;
		jd_tt = jd_utc + ((double)leap_secs + 32.184) / 86400.0;
		jd_ut1 = jd_utc + ut1_utc / 86400.0;
		delta_t = 32.184 + leap_secs - ut1_utc;
		secdif = 0.0;
		tdb2tt(jd_tt, &x, &secdif);
		jd_tdb = jd_tt + secdif / 86400.0;
		jd[0] = jd_tdb;
		jd[1] = 0.0;
		jd0[0] = jd[0] - TA / 86400.0;
		jd0[1] = 0.0;
		//printf("jd[0]：%f jd0[0]：%f\n", jd[0], jd0[0]);
		//printf("theta：%10.9f\n", theta);
		//printf("theta1：%10.9f\n", theta1);
		//printf("thetaE：%f\n", thetaE);
		//printf("thetaM：%f\n",thetaM);
		//printf("deltseta：%f\n", fabs((theta - (thetaE + thetaM))));
	} while (fabs(vpd) > 1e-6);//满月迭代条件
	//} while (OElength < OMlength||theta<theta1);//满月迭代条件
	//} while(OElength < OMlength || fabs((theta - (thetaE + thetaM))) >= 0.00001);//月全食迭代条件
	*/
	//printf("vpd ：%30.29f\n", vpd);
	//jd_utc += 8.0 / 24.0;
	//juliandaytoutc(&jd_utc, &year, &month, &day, &hour1, &min, &sec);
	//printf("月食时间：%f\n", jd_tdb);
	 //printf("日食/月食时间：%d-%d-%d %d:%02d:%02.0f\n", year, month, day, hour1, min, sec);
	//printf("\n\n\n\n\n\n\n \n\n                      2018年农历9月满月时间:%d-%d-%d %d:%02d:%04.2f\n\n\n\n\n\n\n\n \n\n", year, month, day, hour1, min, sec);
    //equ2hor (jd_ut1,delta_t,accuracy,x_pole,y_pole,&obs_loc.on_surf,t_place.ra,t_place.dec,ref_option,&zd,&az,&rar,&decr);
    ephem_close();
   return (0);
}
void juliandaytoutc(double *ju, short int *year, short int *month, short int *day, short int *hour, short int *min, double *second)
{
	int a, b, c, d, e;
	a = (int)(*ju + 0.5);
	b = a + 1537;
	c = (int)((b - 122.1) / 365.25);
	d = (int)(365.25*c);
	e = (int)((b - d) / 30.6001);

	double day1 = b - d - (long)(30.6001*e) + *ju + 0.5 - a;
	*day = (int)day1;
	*month = e - 1 - 12 * (int)(e / 14);
	*year = c - 4715 - (int)((7 + *month) / 10);

	*hour = (int)((day1 - *day)*24.0);
	*min = (int)(((day1 - *day)*24.0 - *hour)*60.0);
	*second = (((day1 - *day) * 24 - *hour) * 60 - *min) * 60;

}
