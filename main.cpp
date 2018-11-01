#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "JTime.h"
#include "MatC.h"

// USER definations
#define USING "KLO0"
#define TRO_CUR "NO"

// STATION definations
#define STA_NAME "UCAL"

// GPS definations
#define F1 1.575420000000000e+09
#define F2 1.227600000000000e+09
#define F5 1.176450000000000e+09
#define F1_F2_2 1.646944444444445
#define F1_F5_2 1.793270321361059

#define L12_M 2.545727780163160
#define L12_N -1.545727780163160
#define L15_M 2.260604327518826
#define L15_N -1.260604327518826
#define L25_M 12.255319148936170
#define L25_N -11.255319148936170

#define GPS_MAX_PRN 32
#define GPS_FLAG 'G'
#define LIGHT_SPEED 299792458.0
#define we 7.2921151467e-5
#define mu 3.986004415E14
#define R1 4.442807633e-10

// RINEX definations
#define END_OF_HEADER "END OF HEADER"
#define MAX_LENGTH_OF_LINE 400
#define SP3_HEADER_LINE 22
#define SP3_EOF "EOF"
#define CLK_AR "AR"
#define CLK_AS "AS"

// solver definations
#define LS_MAX_ITER 10
#define OBS_SIG0    1
#define LS_CONV_THRES 1e-8
#define ACCURACY_CONFIDENCE 100 // meters

struct OBS_FRAME {
	double C1C, L1C, D1C, S1C, C2W, L2W, S2W, C2X, L2X, S2X, C5X, L5X, S5X;
};

struct CLK_FRAME {
	char mode[3];
	char name[5];
	int val_num;
	double CLK_OFF, CLK_DRI;
};

struct BRD_FRAME {
	// PRN / EPOCH / SV CLK
	GPSTime toc;
	double sv_clock_bias;
	double sv_clock_drift;
	double sv_clock_drift_rate;

	// BROADCAST ORBIT - 1
	double idoe_issue_of_data;
	double crs;
	double delta_n;
	double m0;

	// BROADCAST ORBIT - 2
	double cuc;
	double eccentricity;
	double cus;
	double sqrt_a;

	// BROADCAST ORBIT - 3
	double toe;
	double cic;
	double OMEGA;
	double cis;

	// BROADCAST ORBIT - 4
	double i0;
	double crc;
	double omega;
	double omega_dot;

	// BROADCAST ORBIT - 5
	double idot;
	double codes_on_l2;
	double gpsweek;
	double l2_pdata_flag;

	// BROADCAST ORBIT - 6
	double sv_accuracy;
	double sv_health;
	double tgd;
	double iodc_issue_of_data;

	// BROADCAST ORBIT - 7
	double trans_time;
	double fit_interval;
};

// using model Hopefield and NMF
class SimpleTroposphereModel {
	const double Deg = 180.0 / M_PI;
public:
	void set_location(const double * loc)
	{
		location = loc;
	}

	void set_date(int doy) {
		t = cos(2 * M_PI * (doy - 28) / 365.25);
	}

	SimpleTroposphereModel() = default;

	SimpleTroposphereModel(const double * loc, int doy)
	{
		location = loc;
		t = cos(2 * M_PI * (doy - 28) / 365.25);
	}

	double get_dry_std(double elev)
	{
		return get_dry_m(elev) * get_dry_ztd();
	}

	double get_overall_std(double elev)
	{
		return get_dry_std(elev) + get_wet_std(elev);
	}

	double get_wet_std(double elev)
	{
		return get_wet_m(elev) * get_wet_ztd();
	}

	double get_dry_ztd()
	{
		double elev = M_PI / 2;
		double t0, p0, e0, h0;
		double t, p;
		double dend, elev2, hd, rkd;
		double hgt = location[2];

		if (fabs(hgt)>30000.0)   return 0.0;

		t0 = 20 + 273.16;
		p0 = 1013.0;
		e0 = 0.5 * exp(-37.2465 + 0.213166 * t0 - 0.000256908 * t0 * t0);
		h0 = 0;
		t = t0 - 0.0068 * (hgt - h0);
		p = p0 * pow(1.0 - 0.0068 / t0 * (hgt - h0), 5);

		elev2 = elev * elev * Deg * Deg;
		dend = sqrt(elev2 + 6.25) / Deg;

		hd = 148.72 * t0 - 488.3552;
		rkd = 1.552e-5 * p / t * (hd - hgt);
		return  (rkd / sin(dend));
	}

	double get_wet_m(double elev)
	{
		double lat_degree = fabs(location[0] * Deg);
		if (lat_degree > 75)lat_degree = 75;
		else if (lat_degree < 15)lat_degree = 15;
		int lat_index = (int)floor(lat_degree / lat_step) - 1;

		double a, b, c;
		a = wet_abc_ave[0][lat_index] +
			(wet_abc_ave[0][lat_index + 1] - wet_abc_ave[0][lat_index]) *
			(lat_degree - lat_grid[lat_index]) / lat_step;

		b = wet_abc_ave[1][lat_index] +
			(wet_abc_ave[1][lat_index + 1] - wet_abc_ave[1][lat_index]) *
			(lat_degree - lat_grid[lat_index]) / lat_step;

		c = wet_abc_ave[2][lat_index] +
			(wet_abc_ave[2][lat_index + 1] - wet_abc_ave[2][lat_index]) *
			(lat_degree - lat_grid[lat_index]) / lat_step;

		double sinE = sin(elev);
		return (1.0 / (1.0 + a / (1.0 + b / (1.0 + c)))) / (1.0 / (sinE + a / (sinE + b / (sinE + c))));

	}

	double get_wet_ztd()
	{
		double elev = M_PI / 2;
		double t0, p0, e0, h0;
		double t, e;
		double elev2, denw, hw, rkw;
		double hgt = location[2];

		if (fabs(hgt)>30000.0)   return 0.0;

		t0 = 20 + 273.16;
		p0 = 1013.0;
		e0 = 0.5*exp(-37.2465 + 0.213166*t0 - 0.000256908*t0*t0);
		h0 = 0;
		hw = 11000.0;
		t = t0 - 0.0068*(hgt - h0);
		e = e0*pow((1 - 0.0068 / t0*(hgt - h0)), 2.0)*pow((1.0 - (hgt - h0) / hw), 4.0);
		elev2 = elev*elev * Deg * Deg;
		denw = sqrt(elev2 + 2.25) / Deg;

		rkw = 7.46512e-2*(e / t / t)*(hw - hgt);
		return (rkw / sin(denw));
	}

	double get_dry_m(double elev)
	{
		double lat_degree = fabs(location[0] * Deg);
		int lat_index = (int)floor(lat_degree / lat_step) - 1;

		double a, b, c;
		switch (lat_index)
		{
		case -1:
			a = dry_abc_ave[0][0] + dry_abc_ave[0][0] * t;
			b = dry_abc_ave[1][0] + dry_abc_ave[1][0] * t;
			c = dry_abc_ave[2][0] + dry_abc_ave[2][0] * t;
			break;
		case 4:
			a = dry_abc_ave[0][4] + dry_abc_ave[0][4] * t;
			b = dry_abc_ave[1][4] + dry_abc_ave[1][4] * t;
			c = dry_abc_ave[2][4] + dry_abc_ave[2][4] * t;
			break;
		default:
			a = dry_abc_ave[0][lat_index] +
				(dry_abc_ave[0][lat_index + 1] - dry_abc_ave[0][lat_index]) *
				(lat_degree - lat_grid[lat_index]) / lat_step;
			a += dry_abc_amp[0][lat_index] +
				(dry_abc_amp[0][lat_index + 1] - dry_abc_amp[0][lat_index]) *
				(lat_degree - lat_grid[lat_index]) / lat_step * t;

			b = dry_abc_ave[1][lat_index] +
				(dry_abc_ave[1][lat_index + 1] - dry_abc_ave[1][lat_index]) *
				(lat_degree - lat_grid[lat_index]) / lat_step;
			b += dry_abc_amp[1][lat_index] +
				(dry_abc_amp[1][lat_index + 1] - dry_abc_amp[1][lat_index]) *
				(lat_degree - lat_grid[lat_index]) / lat_step * t;

			c = dry_abc_ave[2][lat_index] +
				(dry_abc_ave[2][lat_index + 1] - dry_abc_ave[2][lat_index]) *
				(lat_degree - lat_grid[lat_index]) / lat_step;
			c += dry_abc_amp[2][lat_index] +
				(dry_abc_amp[2][lat_index + 1] - dry_abc_amp[2][lat_index]) *
				(lat_degree - lat_grid[lat_index]) / lat_step * t;
			break;

		}
		double sinE = sin(elev);
		return (1.0 / (1.0 + a / (1.0 + b / (1 + c)))) / (1.0 / (sinE + a / (sinE + b / (sinE + c)))) +
			(1.0 / sinE - (1.0 / (1.0 + abc_ht[0] / (1.0 + abc_ht[1] / (1.0 + abc_ht[2])))) / (1.0 / (sinE + abc_ht[0] / (sinE + abc_ht[1] / (sinE + abc_ht[2]))))) * location[2] / 1000;
	}
protected:
	const double * location;
	double t;


	const double abc_ht[3]{ 2.53E-5, 5.49E-3,1.14E-3 };
	const int lat_step = 15;
	const int lat_grid[5]{ 15, 30,45, 60,75 };

	const double dry_abc_ave[3][5]{
		{ 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3 },
		{ 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3 },
		{ 62.620505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3 }
	};

	const double dry_abc_amp[3][5]{
		{ 0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5 },
		{ 0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5 },
		{ 0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5 }
	};

	const double wet_abc_ave[3][5]{
		{ 5.8021879E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4 },
		{ 1.4275268E-3, 1.5138625E-3, 1.4572572E-3, 1.5007428E-3, 1.7599082E-3 },
		{ 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736039E-2 }
	};
private:

};

struct SP3_FRAME {
	double X, Y, Z, dt;
};
const char * brd_files[]{
	{ "./RAW/brdc0080.18n" },
	{ "./RAW/brdc0090.18n" },
	{ "./RAW/brdc0100.18n" },
	{ "./RAW/brdc0110.18n" },
	{ "./RAW/brdc0120.18n" },
	{ "./RAW/brdc0130.18n" },
	{ "./RAW/brdc0140.18n" },
	{ "./RAW/brdc0150.18n" },
	{ "./RAW/brdc0160.18n" },
	{ "./RAW/brdc0170.18n" },
	{ "./RAW/brdc0180.18n" },
	{ "./RAW/brdc0190.18n" },
	{ "./RAW/brdc0200.18n" },
	{ "./RAW/brdc0210.18n" },
};

const char * obs_files[]{
	{ "./RAW/UCAL00CAN_S_20180080000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180090000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180100000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180110000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180120000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180130000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180140000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180150000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180160000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180170000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180180000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180190000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180200000_01D_30S_MO.rnx" },
	{ "./RAW/UCAL00CAN_S_20180210000_01D_30S_MO.rnx" }
};
const char * sp3_files[]{
	{ "./RAW/igs19831.sp3" },
	{ "./RAW/igs19832.sp3" },
	{ "./RAW/igs19833.sp3" },
	{ "./RAW/igs19834.sp3" },
	{ "./RAW/igs19835.sp3" },
	{ "./RAW/igs19836.sp3" },
	{ "./RAW/igs19840.sp3" },
	{ "./RAW/igs19841.sp3" },
	{ "./RAW/igs19842.sp3" },
	{ "./RAW/igs19843.sp3" },
	{ "./RAW/igs19844.sp3" },
	{ "./RAW/igs19845.sp3" },
	{ "./RAW/igs19846.sp3" },
	{ "./RAW/igs19850.sp3" },
};
const char * clk_files[]{
	{ "./RAW/igs19831.clk_30s" },
	{ "./RAW/igs19832.clk_30s" },
	{ "./RAW/igs19833.clk_30s" },
	{ "./RAW/igs19834.clk_30s" },
	{ "./RAW/igs19835.clk_30s" },
	{ "./RAW/igs19836.clk_30s" },
	{ "./RAW/igs19840.clk_30s" },
	{ "./RAW/igs19841.clk_30s" },
	{ "./RAW/igs19842.clk_30s" },
	{ "./RAW/igs19843.clk_30s" },
	{ "./RAW/igs19844.clk_30s" },
	{ "./RAW/igs19845.clk_30s" },
	{ "./RAW/igs19846.clk_30s" },
	{ "./RAW/igs19850.clk_30s" },
};
// consts level
const double station_xyz[3] = { -1641945.2000, -3664804.1000, 4940009.3000 };
const double station_blh[3] = { 0.891513720275876, -1.99201147820819, 1118.80358502921 };
const double center_earth[3] = { 0,0,0 };
const int num_of_date = 14;
double sinB = sin(station_blh[0]);
double cosB = cos(station_blh[0]);
double sinL = sin(station_blh[1]);
double cosL = cos(station_blh[1]);

// file input level
UTC utc_obs = UTC::current_utc();
UTC utc_sp3 = UTC::current_utc();
UTC utc_clk = UTC::current_utc();
UTC utc_brd = UTC::current_utc();

OBS_FRAME OBS[GPS_MAX_PRN];
BRD_FRAME BRD[GPS_MAX_PRN];
SP3_FRAME BRP[GPS_MAX_PRN];
SP3_FRAME SP3[GPS_MAX_PRN];
CLK_FRAME CLK[GPS_MAX_PRN];
CLK_FRAME REC;

bool obs_available[GPS_MAX_PRN] = { false };
bool sp3_available[GPS_MAX_PRN] = { false };
bool clk_available[GPS_MAX_PRN] = { false };
bool brd_available[GPS_MAX_PRN] = { false };
bool rec_available = false;

int sat_num = 0;

FILE * ofp, *sfp, *cfp, *bfp;
char line_buffer[MAX_LENGTH_OF_LINE] = "";
char nav_buffer[20] = "";

char SVN[4] = "";
CLK_FRAME temp_clk;
double seconds = 0.0;
int prn_check = 0;
int brd_prn = 0;

// solver level
SimpleTroposphereModel tro_model;
OBS_FRAME * S_OBS[GPS_MAX_PRN];
SP3_FRAME * S_SP3[GPS_MAX_PRN];
BRD_FRAME * S_BRD[GPS_MAX_PRN];
SP3_FRAME * S_BRP[GPS_MAX_PRN];
CLK_FRAME * S_CLK[GPS_MAX_PRN];

int PRN[GPS_MAX_PRN];
bool solve_available[GPS_MAX_PRN];
double solution[4];
double correction_values[GPS_MAX_PRN];
double ion_corrections[GPS_MAX_PRN];

// analysis level
GPSTime gpst;
double sat_elevation[GPS_MAX_PRN];
double ref_distance[GPS_MAX_PRN];
double sat_azimuth[GPS_MAX_PRN];
double sat_ek[GPS_MAX_PRN];
double ion_param[8];
double klo[3][GPS_MAX_PRN];
double com[3][GPS_MAX_PRN];

FILE * outfp = fopen("out.txt", "w");
FILE * klo1fp = fopen("klo1.txt", "w");
FILE * klo2fp = fopen("klo2.txt", "w");
FILE * klo5fp = fopen("klo5.txt", "w");
FILE * com12fp = fopen("com12.txt", "w");
FILE * com15fp = fopen("com15.txt", "w");
FILE * com25fp = fopen("com25.txt", "w");

double klobuchar(int index)
{
	double I = 0.0;
	double sai = 0.0137 / (sat_elevation[index] / M_PI + 0.11) - 0.022;
	double phyi = station_blh[0] / M_PI + sai * cos(sat_azimuth[index]);
	if (phyi > 0.416) phyi = 0.416;
	if (phyi < -0.416) phyi = -0.416;
	double namdai = station_blh[1] / M_PI + sai * sin(sat_azimuth[index]) / cos(phyi * M_PI);
	double phym = phyi + 0.064 * cos((namdai - 1.617) * M_PI);
	double t = 43200.0 * namdai + gpst.sec;
	if (t >= 86400) t = fmod(t, 86400.0);
	else if (t < 0) t = t + 86400;

	double Ai = ion_param[0] + ion_param[1] * phym + ion_param[2] * phym * phym + ion_param[3] * phym * phym * phym;
	if (Ai < 0) Ai = 0.0;
	double Pi = ion_param[4] + ion_param[5] * phym + ion_param[6] * phym * phym + ion_param[7] * phym * phym * phym;
	if (Pi < 72000) Pi = 72000;

	double Xi = 2 * M_PI * (t - 50400) / Pi;
	double f = 1.0 + 16.0 * (0.53 - sat_elevation[index] / M_PI) * (0.53 - sat_elevation[index] / M_PI) * (0.53 - sat_elevation[index] / M_PI);
	double X = 1.0 - Xi * Xi / 2.0 + Xi * Xi * Xi * Xi / 24.0;
	if (fabs(Xi) <= 1.57) {
		I = (5.0e-9 + Ai * X) * f;
	}
	else
		I = 5.0e-9 * f;
	return I;
}

double distance(const double * p1, const double * p2, int dim = 3)
{
	double tot = 0;
	for (int i = 0; i < dim; i++)
	{
		tot += (p1[i] - p2[i]) * (p1[i] - p2[i]);
	}
	return sqrt(tot);
}

void elevation_and_azimuth(SP3_FRAME * satellite, int index)
{
	double dpos[3] = { 0 };
	double ori[3]{ 0,0,0 };
	dpos[0] = satellite->X - station_xyz[0];
	dpos[1] = satellite->Y - station_xyz[1];
	dpos[2] = satellite->Z - station_xyz[2];

	double user_distance_to_earth = distance(station_xyz, center_earth);
	double mod = sqrt(dpos[0] * dpos[0] + dpos[1] * dpos[1] + dpos[2] * dpos[2]);
	if (fabs(user_distance_to_earth * mod < 1.0)) {
		sat_elevation[index] = M_PI_2;
	}
	else {
		double m = dpos[0] * station_xyz[0] + dpos[1] * station_xyz[1] + dpos[2] * station_xyz[2];
		double n = m / (mod * user_distance_to_earth);
		sat_elevation[index] = M_PI_2 - acos(n);
	}

	double N = -sinB * cosL * dpos[0] - sinB * sinL * dpos[1] + cosB * dpos[2];
	double E = -sinL * dpos[0] + cosL* dpos[1];
	sat_azimuth[index] = atan2(E, N);
}


void fetch_clk(FILE * fp)
{
	while (!utc_clk.equals(&utc_sp3)) {
		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		sscanf(line_buffer, "%s%s%d%d%d%d%d%lf%d%lf%lf",
			&temp_clk.mode, &temp_clk.name,
			&utc_clk.year, &utc_clk.month, &utc_clk.date, &utc_clk.hour, &utc_clk.minute, &seconds,
			&temp_clk.val_num, &temp_clk.CLK_OFF, &temp_clk.CLK_DRI
		);
		utc_clk.year = utc_clk.year - 2000;
		utc_clk.sec = (int)round(seconds);
	}

	//long long origin = ftell(fp);
	for (;; ) {
		
		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		sscanf(line_buffer, "%s%s%d%d%d%d%d%lf%d%lf%lf",
			&temp_clk.mode, &temp_clk.name,
			&utc_clk.year, &utc_clk.month, &utc_clk.date, &utc_clk.hour, &utc_clk.minute, &seconds,
			&temp_clk.val_num, &temp_clk.CLK_OFF, &temp_clk.CLK_DRI
		);
		utc_clk.year = utc_clk.year - 2000;
		if (seconds == utc_clk.sec)
		{
			if (strcmp(temp_clk.mode, CLK_AR) == 0) {
				if (strcmp(temp_clk.name, STA_NAME) == 0) {
					memcpy(&REC, &temp_clk, sizeof(CLK_FRAME));
					rec_available = true;
				}
			}
			else if(strcmp(temp_clk.mode, CLK_AS) == 0){
				if (sscanf(temp_clk.name, "G%d", &prn_check)) {
					memcpy(CLK + prn_check - 1, &temp_clk, sizeof(CLK_FRAME));
					clk_available[prn_check - 1] = true;
				}
			}
		}
		else break;
	}

	fseek(fp, -20 - strlen(line_buffer), SEEK_CUR);
}

void fetch_sp3(FILE * fp)
{
	fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
	if (strncmp(line_buffer, SP3_EOF, 3) == 0) throw -1;

	sscanf(line_buffer, "*  %d %d %d %d %d %lf",
		&utc_sp3.year, &utc_sp3.month, &utc_sp3.date, &utc_sp3.hour, &utc_sp3.minute, &seconds);
	utc_sp3.year = utc_sp3.year - 2000;
	utc_sp3.sec = (int)round(seconds);

	for (int i = 0; i < GPS_MAX_PRN; i++)
	{
		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		double * arr = (double*)(SP3 + i);
		sscanf(line_buffer, "PG%d %lf %lf %lf %lf", 
			&prn_check, arr, arr + 1, arr + 2, arr + 3);
		for (int j = 0; j < 3; j++)arr[j] *= 1e3;

		sp3_available[i] = true;
	}
}

inline double _fastcall extract_double(const char * pointee)
{
	if (strlen(pointee) < 19) return 0;
	strncpy(nav_buffer, pointee, 19);
	nav_buffer[15] = 'E';
	return atof(nav_buffer);
}

void brd_grab_time(FILE * fp)
{
	if (feof(fp)) throw -1;
	double sec;
	fscanf(fp, "%d%d%d%d%d%d%lf",
		&brd_prn,
		&utc_brd.year,
		&utc_brd.month,
		&utc_brd.date,
		&utc_brd.hour,
		&utc_brd.minute,
		&sec
	);
	utc_brd.sec = (int)round(sec);
}

void fetch_brd(FILE * fp)
{
	while (utc_obs.larger_than(utc_brd))
	{
		brd_prn--;
		// PRN / EPOCH / SV CLK
		if (fgets(line_buffer, MAX_LENGTH_OF_LINE, fp) == NULL) throw -1;
		BRD[brd_prn].toc = GPSTime(utc_brd);
		BRD[brd_prn].sv_clock_bias = extract_double(line_buffer);
		BRD[brd_prn].sv_clock_drift = extract_double(line_buffer + 19);
		BRD[brd_prn].sv_clock_drift_rate = extract_double(line_buffer + 38);

		// BROADCAST ORBIT - 1
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD[brd_prn].idoe_issue_of_data = extract_double(line_buffer);
		BRD[brd_prn].crs = extract_double(line_buffer + 19);
		BRD[brd_prn].delta_n = extract_double(line_buffer + 38);
		BRD[brd_prn].m0 = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 2
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD[brd_prn].cuc = extract_double(line_buffer);
		BRD[brd_prn].eccentricity = extract_double(line_buffer + 19);
		BRD[brd_prn].cus = extract_double(line_buffer + 38);
		BRD[brd_prn].sqrt_a = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 3
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD[brd_prn].toe = extract_double(line_buffer);
		BRD[brd_prn].cic = extract_double(line_buffer + 19);
		BRD[brd_prn].OMEGA = extract_double(line_buffer + 38);
		BRD[brd_prn].cis = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 4
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD[brd_prn].i0 = extract_double(line_buffer);
		BRD[brd_prn].crc = extract_double(line_buffer + 19);
		BRD[brd_prn].omega = extract_double(line_buffer + 38);
		BRD[brd_prn].omega_dot = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 5
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD[brd_prn].idot = extract_double(line_buffer);
		BRD[brd_prn].codes_on_l2 = extract_double(line_buffer + 19);
		BRD[brd_prn].gpsweek = extract_double(line_buffer + 38);
		BRD[brd_prn].l2_pdata_flag = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 6
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD[brd_prn].sv_accuracy = extract_double(line_buffer);
		BRD[brd_prn].sv_health = extract_double(line_buffer + 19);
		BRD[brd_prn].tgd = extract_double(line_buffer + 38);
		BRD[brd_prn].iodc_issue_of_data = extract_double(line_buffer + 57);

		// BROADCAST ORBIT - 7
		fgets(line_buffer, 4, fp);
		fgets(line_buffer, 90, fp);
		BRD[brd_prn].trans_time = extract_double(line_buffer);
		BRD[brd_prn].fit_interval = extract_double(line_buffer + 19);

		brd_available[brd_prn] = true;

		brd_grab_time(fp);
	}
}

inline double _fastcall get_atan(double z, double y)
{
	double x = 0;
	if (z == 0)x = M_PI / 2;
	else if (y == 0)x = M_PI;
	else {
		x = atan(fabs(y / z));
		if ((y > 0) && (z < 0))x = M_PI - x;
		else if ((y < 0) && (z < 0))x = M_PI + x;
		else if ((y < 0) && (z > 0))x = 2 * M_PI - x;
	}
	return x;
}

bool brdc_satell(int index)
{
	BRD_FRAME * f = S_BRD[index];
	double n0 = sqrt(mu) / pow(f->sqrt_a, 3);
	double n = n0 + f->delta_n;

	double t = gpst.sec - S_OBS[index]->C1C / LIGHT_SPEED;
	double tk = t - f->toe;
	if (fabs(tk)>7200)
		return false;
	if (tk > 302400)
		tk -= 604800;
	else if (tk < -302400)
		tk += 604800;

	double Mk = f->m0 + n * tk;

	sat_ek[index] = Mk;
	double Ek2 = 0;
	while (1)
	{
		Ek2 = Mk + f->eccentricity * sin(sat_ek[index]);
		if (fabs(sat_ek[index] - Ek2) <= 1.0e-12)break;
		sat_ek[index] = Ek2;
	}

	double sqt_1_e2 = sqrt(1 - pow(f->eccentricity, 2));
	double cosfk = (cos(sat_ek[index]) - f->eccentricity) / (1 - f->eccentricity * cos(sat_ek[index]));
	double sinfk = (sqt_1_e2 * sin(sat_ek[index])) / (1 - f->eccentricity * cos(sat_ek[index]));
	double fk = get_atan((cos(sat_ek[index]) - f->eccentricity), sqt_1_e2 * sin(sat_ek[index]));
	
	double faik = fk + f->omega;

	double cos_2_faik = cos(2 * faik);
	double sin_2_faik = sin(2 * faik);
	double su = f->cuc * cos_2_faik + f->cus * sin_2_faik;
	double sr = f->crc * cos_2_faik + f->crs * sin_2_faik;
	double si = f->cic * cos_2_faik + f->cis * sin_2_faik;

	double uk = faik + su;
	double rk = f->sqrt_a * f->sqrt_a * (1 - f->eccentricity * cos(sat_ek[index])) + sr;
	double ik = f->i0 + si + f->idot * tk;

	double xk = rk * cos(uk);
	double yk = rk * sin(uk);
	double zk = 0;

	double L = f->OMEGA + (f->omega_dot - we) * tk - we * f->toe;

	S_BRP[index]->X = xk * cos(L) - yk * cos(ik) * sin(L);
	S_BRP[index]->Y = xk * sin(L) + yk * cos(ik) * cos(L);
	S_BRP[index]->Z = yk * sin(ik);

	return true;
}

void fetch_obs(FILE * fp)
{
	// time & sat_num
	int health = 0;

	fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
	sscanf(line_buffer, "> %d %d %d %d %d %lf %d %d",
		&utc_obs.year, &utc_obs.month, &utc_obs.date, &utc_obs.hour, &utc_obs.minute, &seconds, &health, &sat_num);
	utc_obs.year = utc_obs.year - 2000;
	utc_obs.sec = (int)round(seconds);

	//while (!utc_obs.equals(&utc_sp3))
	//{
	//	for (int i = 0; i < sat_num; i++)
	//	{
	//		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
	//	}

	//	fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
	//	sscanf(line_buffer, "> %d %d %d %d %d %lf %d %d",
	//		&utc_obs.year, &utc_obs.month, &utc_obs.date, &utc_obs.hour, &utc_obs.minute, &seconds, &health, &sat_num);
	//	utc_obs.year = utc_obs.year - 2000;
	//	utc_obs.sec = (int)round(seconds);
	//}

	for (int i = 0; i < sat_num; i++)
	{
		fgets(SVN, 4, fp);
		if (SVN[0] == GPS_FLAG)
		{
			int prn = atoi(SVN + 1);
			int index = 0;
			double * ptr = (double*)(OBS + prn - 1);
			while (fscanf(fp, "%lf", ptr + index)) index++;
			obs_available[prn - 1] = true;
		}
		else {
			fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		}
	}
}

bool skip_obs_header(FILE * fp)
{
	while (!feof(fp))
	{
		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		if (strncmp(line_buffer + 60, END_OF_HEADER, 13) == 0)
			return true;
	}
	return false;
}

bool skip_brd_header(FILE * fp)
{
	for (int i = 0; i < 4; i++)
		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
	
	for (int i = 0; i < 60; i++)
		if (line_buffer[i] == 'D')line_buffer[i] = 'E';
	
	sscanf(line_buffer, "%lf%lf%lf%lf",
		ion_param, ion_param + 1, ion_param + 2, ion_param + 3);

	fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);

	for (int i = 0; i < 60; i++)
		if (line_buffer[i] == 'D')line_buffer[i] = 'E';

	sscanf(line_buffer, "%lf%lf%lf%lf",
		ion_param + 4, ion_param + 5, ion_param + 6, ion_param + 7);

	if (skip_obs_header(fp))
		brd_grab_time(fp);
	else return false;
	return true;
}

bool skip_clk_header(FILE * fp)
{
	while (!feof(fp))
	{
		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		if (strncmp(line_buffer + 60, END_OF_HEADER, 13) == 0)
			return true;
	}
	return false;
}

bool skip_sp3_header(FILE * fp)
{
	for (int i = 0; i < SP3_HEADER_LINE; i++)
	{
		if (feof(fp)) return false;
		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
	}
	return true;
}

inline double _fastcall distance(SP3_FRAME * f1, SP3_FRAME * f2)
{
	return sqrt((f1->X - f2->X)*(f1->X - f2->X) + (f1->Y - f2->Y) * (f1->Y - f2->Y) + (f1->Z - f2->Z) * (f1->Z - f2->Z));
}

void corrections(int index)
{
	double dtc = gpst.minus(&S_BRD[index]->toc) - S_OBS[index]->C1C / LIGHT_SPEED;
	double s = S_BRD[index]->sv_clock_bias
		+ S_BRD[index]->sv_clock_drift      * dtc
		+ S_BRD[index]->sv_clock_drift_rate * dtc;

	double r = (R1 * S_BRD[index]->eccentricity * sin(sat_ek[index]) * S_BRD[index]->sqrt_a);

	
	correction_values[PRN[index] - 1] -= r * LIGHT_SPEED;                  //相对论
	correction_values[PRN[index] - 1] -= S_BRD[index]->tgd * LIGHT_SPEED;  //群延迟

	if (strcmp(TRO_CUR, "YES") == 0) // 对流层
		correction_values[PRN[index] - 1] -= tro_model.get_overall_std(sat_elevation[PRN[index] - 1]);

	if (S_OBS[index]->C1C)S_OBS[index]->C1C += correction_values[PRN[index] - 1];
	if (S_OBS[index]->C2X)S_OBS[index]->C2X += correction_values[PRN[index] - 1];
	if (S_OBS[index]->C5X)S_OBS[index]->C5X += correction_values[PRN[index] - 1];

	double dA = we * (S_OBS[sat_num]->C1C / LIGHT_SPEED);
	double Xc = S_BRP[sat_num]->Y * dA;
	double Yc = S_BRP[sat_num]->X * dA;

	S_BRP[sat_num]->X += Xc;
	S_BRP[sat_num]->Y -= Yc;
}

bool valid_check(int i)
{
	if (strcmp(USING, "COM0") == 0 && (OBS[i].C1C == 0 || OBS[i].C2X == 0))
		return false;
	else if (strcmp(USING, "COM1") == 0 && (OBS[i].C1C == 0 || OBS[i].C5X == 0))
		return false;
	else if (strcmp(USING, "COM2") == 0 && (OBS[i].C2X == 0 || OBS[i].C5X == 0))
		return false;

	else if (strcmp(USING, "KLO0") == 0 && OBS[i].C1C == 0)
		return false;
	else if (strcmp(USING, "KLO1") == 0 && OBS[i].C2X == 0)
		return false;
	else if (strcmp(USING, "KLO2") == 0 && OBS[i].C5X == 0)
		return false;

	else if (strcmp(USING, "C1C") == 0 && OBS[i].C1C == 0)
		return false;
	else if (strcmp(USING, "C2X") == 0 && OBS[i].C2X == 0)
		return false;
	else if (strcmp(USING, "C5X") == 0 && OBS[i].C5X == 0)
		return false;
	else if (strcmp(USING, "C2W") == 0 && OBS[i].C2W == 0)
		return false;

	return true;
}

void overall_check()
{
	sat_num = 0;
	gpst = GPSTime(utc_obs);

	for (int i = 0; i < GPS_MAX_PRN; i++)
	{
		if (obs_available[i] && brd_available[i])
		{
			if (valid_check(i))
			{
				S_OBS[sat_num] = OBS + i;
				S_BRD[sat_num] = BRD + i;
				S_BRP[sat_num] = BRP + i;
				PRN[sat_num] = i + 1;

				brdc_satell(sat_num);
				elevation_and_azimuth(SP3 + i, i);
				corrections(sat_num);
				
				solve_available[i] = true;
				sat_num++;
			}
		}
	}
}

void erace()
{
	memset(klo[0], 0, sizeof(double) * GPS_MAX_PRN);
	memset(klo[1], 0, sizeof(double) * GPS_MAX_PRN);
	memset(klo[2], 0, sizeof(double) * GPS_MAX_PRN);
	memset(com[0], 0, sizeof(double) * GPS_MAX_PRN);
	memset(com[1], 0, sizeof(double) * GPS_MAX_PRN);
	memset(com[2], 0, sizeof(double) * GPS_MAX_PRN);

	memset(correction_values, 0, sizeof(double) * GPS_MAX_PRN);
	memset(ion_corrections, 0, sizeof(double) * GPS_MAX_PRN);

	memset(OBS, 0, sizeof(OBS_FRAME) * GPS_MAX_PRN);
	memset(SP3, 0, sizeof(SP3_FRAME) * GPS_MAX_PRN);
	memset(CLK, 0, sizeof(CLK_FRAME) * GPS_MAX_PRN);
	memset(BRP, 0, sizeof(CLK_FRAME) * GPS_MAX_PRN);

	memset(S_OBS, 0, sizeof(OBS_FRAME*) * GPS_MAX_PRN);
	memset(S_SP3, 0, sizeof(SP3_FRAME*) * GPS_MAX_PRN);
	memset(S_CLK, 0, sizeof(CLK_FRAME*) * GPS_MAX_PRN);
	memset(S_BRP, 0, sizeof(CLK_FRAME*) * GPS_MAX_PRN);
	memset(S_BRD, 0, sizeof(BRD_FRAME*) * GPS_MAX_PRN);

	memset(obs_available, 0, sizeof(bool) * GPS_MAX_PRN);
	memset(sp3_available, 0, sizeof(bool) * GPS_MAX_PRN);
	memset(clk_available, 0, sizeof(bool) * GPS_MAX_PRN);
	memset(solve_available, 0, sizeof(bool) * GPS_MAX_PRN);

	sat_num = 0;
	memset(solution, 0, sizeof(double) * 4);

	memset(sat_elevation, 0, sizeof(double) * GPS_MAX_PRN);
	memset(sat_azimuth, 0, sizeof(double) * GPS_MAX_PRN);
	memset(ref_distance, 0, sizeof(double) * GPS_MAX_PRN);

	memset(line_buffer, 0, sizeof(char) * MAX_LENGTH_OF_LINE);
	memset(SVN, 0, sizeof(char) * 4);
	seconds = 0.0;
	prn_check = 0;
}
double rough_dis = 0.0;

bool rough_solve()
{
	Matrix * L = malloc_mat(sat_num, 1);
	Matrix * A = malloc_mat(sat_num, 4);
	Matrix * Cl = malloc_mat(sat_num, sat_num); 
	Matrix * r = malloc_mat(sat_num, 1);
	Matrix * Q = malloc_mat(4, 4);
	Matrix * δ = malloc_mat(4, 1);
	Matrix * Cx = malloc_mat(4, 4);

	double * DX0 = (double*)alloca(sat_num * sizeof(double)); 
	double * DY0 = (double*)alloca(sat_num * sizeof(double));
	double * DZ0 = (double*)alloca(sat_num * sizeof(double));
	double * S = (double*)alloca(sat_num * sizeof(double)); 

	double last_solution[4] = { 0,0,0,0 };
	memcpy(solution, station_xyz, sizeof(double) * 3);

	for (int i = 0; i < LS_MAX_ITER; i++)
	{
		memcpy(last_solution, solution, sizeof(double) * 4);

		for (int j = 0; j < sat_num; j++)
		{
			DX0[j] = S_BRP[j]->X - solution[0];
			DY0[j] = S_BRP[j]->Y - solution[1];
			DZ0[j] = S_BRP[j]->Z - solution[2];
			S[j] = sqrt(DX0[j] * DX0[j] + DY0[j] * DY0[j] + DZ0[j] * DZ0[j]);
		}

		// get A, L, Cl matrices.
		for (int j = 0; j < sat_num; j++)
		{
			double sinval = sin(sat_elevation[PRN[j] - 1]);
			Cl->data[j][j] = OBS_SIG0 * OBS_SIG0 / (sinval * sinval);

			// for common spp
			L->data[j][0] = S_OBS[j]->C1C - S[j] - solution[3] + S_BRD[j]->sv_clock_bias * LIGHT_SPEED;

			// A matrix
			A->data[j][0] = -DX0[j] / S[j];  // for X
			A->data[j][1] = -DY0[j] / S[j];  // for Y
			A->data[j][2] = -DZ0[j] / S[j];  // for Z
			A->data[j][3] = 1;               // for cdt
		}

		LMS(L, A, Cl, δ, Q, r, Cx);   

		for (int j = 0; j < 4; j++)
			solution[j] += δ->data[j][0];

		// If converged
		if (distance(last_solution, solution, 3) <= LS_CONV_THRES) {
			rough_dis = distance(solution, station_xyz, 3);
			if(rough_dis <= ACCURACY_CONFIDENCE)
				return true;
			return false;
		}
	}

	return false;
}

void klo_calculation()
{
	for (int i = 0; i < sat_num; i++)
	{
		double k_f1 = klobuchar(PRN[i] - 1) * LIGHT_SPEED;
		if (S_OBS[i]->C1C) klo[0][PRN[i] - 1] = S_OBS[i]->C1C - k_f1;
		if (S_OBS[i]->C2X) klo[1][PRN[i] - 1] = S_OBS[i]->C2X - k_f1 * F1_F2_2;
		if (S_OBS[i]->C5X) klo[2][PRN[i] - 1] = S_OBS[i]->C5X - k_f1 * F1_F5_2;
	}
}

void com_calculation()
{
	for (int i = 0; i < sat_num; i++)
	{
		if (S_OBS[i]->C1C && S_OBS[i]->C2X) 
			com[0][PRN[i] - 1] = S_OBS[i]->C1C * L12_M + S_OBS[i]->C2X * L12_N;
		if (S_OBS[i]->C1C && S_OBS[i]->C5X) 
			com[1][PRN[i] - 1] = S_OBS[i]->C1C * L15_M + S_OBS[i]->C5X * L15_N;
		if (S_OBS[i]->C2X && S_OBS[i]->C5X) 
			com[2][PRN[i] - 1] = S_OBS[i]->C2X * L25_M + S_OBS[i]->C5X * L25_N;
	}
}

void output_solution(FILE * fp)
{
	fprintf(fp, "%d\t%lf\t%lf\t%lf\t%lf\n",
		gpst.sec, solution[0], solution[1], solution[2], solution[3]);
}

void output_correction()
{
	for (int i = 0; i < GPS_MAX_PRN; i++)
	{
		fprintf(klo1fp, "%lf\t", klo[0][i]);
		fprintf(klo2fp, "%lf\t", klo[1][i]);
		fprintf(klo5fp, "%lf\t", klo[2][i]);

		fprintf(com12fp, "%lf\t", com[0][i]);
		fprintf(com15fp, "%lf\t", com[1][i]);
		fprintf(com25fp, "%lf\t", com[2][i]);
	}
}

bool make_mode(const char* mode)
{
	if (strcmp(mode, "COM0") == 0)
		for (int i = 0; i < sat_num; i++)
			if (com[0][PRN[i] - 1])
				ion_corrections[PRN[i] - 1] = S_OBS[i]->C1C - com[0][PRN[i] - 1], 
				S_OBS[i]->C1C = com[0][PRN[i] - 1];
			else return false;
	else if (strcmp(mode, "COM1") == 0)
		for (int i = 0; i < sat_num; i++)
			if (com[1][PRN[i] - 1])
				ion_corrections[PRN[i] - 1] = S_OBS[i]->C1C - com[1][PRN[i] - 1],
				S_OBS[i]->C1C = com[1][PRN[i] - 1];
			else return false;
	else if (strcmp(mode, "COM2") == 0)
		for (int i = 0; i < sat_num; i++)
			if (com[2][PRN[i] - 1])
				ion_corrections[PRN[i] - 1] = S_OBS[i]->C1C - com[2][PRN[i] - 1],
				S_OBS[i]->C1C = com[2][PRN[i] - 1];
			else return false;
	else if (strcmp(mode, "KLO0") == 0)
		for (int i = 0; i < sat_num; i++)
			if (klo[0][PRN[i] - 1])
				ion_corrections[PRN[i] - 1] = S_OBS[i]->C1C - klo[0][PRN[i] - 1],
				S_OBS[i]->C1C = klo[0][PRN[i] - 1];
			else return false;
	else if (strcmp(mode, "KLO1") == 0)
		for (int i = 0; i < sat_num; i++)
			if (klo[1][PRN[i] - 1])
				ion_corrections[PRN[i] - 1] = S_OBS[i]->C1C - klo[1][PRN[i] - 1],
				S_OBS[i]->C1C = klo[1][PRN[i] - 1];
			else return false;
	else if (strcmp(mode, "KLO2") == 0)
		for (int i = 0; i < sat_num; i++)
			if (klo[2][PRN[i] - 1])
				ion_corrections[PRN[i] - 1] = S_OBS[i]->C1C - klo[2][PRN[i] - 1],
				S_OBS[i]->C1C = klo[2][PRN[i] - 1];
			else return false;
	else if (strcmp(mode, "C2X") == 0)
		for (int i = 0; i < sat_num; i++)
			if (S_OBS[i]->C2X)
				ion_corrections[PRN[i] - 1] = S_OBS[i]->C1C - S_OBS[i]->C2X,
				S_OBS[i]->C1C = S_OBS[i]->C2X;
			else return false;
	else if (strcmp(mode, "C2W") == 0)
		for (int i = 0; i < sat_num; i++)
			if (S_OBS[i]->C2W)
				ion_corrections[PRN[i] - 1] = S_OBS[i]->C1C - S_OBS[i]->C2W,
				S_OBS[i]->C1C = S_OBS[i]->C2W;
			else return false;
	else if (strcmp(mode, "C5X") == 0)
		for (int i = 0; i < sat_num; i++)
			if (S_OBS[i]->C5X)
				ion_corrections[PRN[i] - 1] = S_OBS[i]->C1C - S_OBS[i]->C5X,
				S_OBS[i]->C1C = S_OBS[i]->C5X;
			else return false;
	
	return sat_num > 3;
}

void execute(int index)
{
	ofp = fopen(obs_files[index], "r");
	sfp = fopen(sp3_files[index], "r");
	cfp = fopen(clk_files[index], "r");
	bfp = fopen(brd_files[index], "r");

	skip_obs_header(ofp);
	skip_sp3_header(sfp);
	skip_clk_header(cfp);
	skip_brd_header(bfp);

	//tro_model.set_location(station_blh);
	
	while (!feof(ofp) && !feof(sfp) && !feof(cfp) && !feof(bfp))
	{

		try {
			erace();
			//fetch_sp3(sfp);
			//fetch_clk(cfp);
			fetch_obs(ofp);
			fetch_brd(bfp);
			//tro_model.set_date(utc_obs.get_doy());

			overall_check();
			
			klo_calculation();
			com_calculation();

			if (!make_mode(USING)) continue;

			if (rough_solve())
			{
				output_solution(outfp);
				output_correction();
			}
		}
		catch (...)
		{
			break;
		}
	}

}

int main()
{
	for (int i = 0; i < num_of_date; i++)
	{
		execute(i);
	}
	_fcloseall();
}

