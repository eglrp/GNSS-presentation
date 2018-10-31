
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "JTime.h"

// STATION definations
#define STA_NAME "UCAL"

// GPS definations
#define GPS_MAX_PRN 32
#define GPS_FLAG 'G'

// RINEX definations
#define END_OF_HEADER "END OF HEADER"
#define MAX_LENGTH_OF_LINE 400
#define SP3_HEADER_LINE 22
#define CLK_AR "AR"
#define CLK_AS "AS"

struct OBS_FRAME {
	double C1C, L1C, D1C, S1C, C2W, L2W, S2W, C2X, L2X, S2X, C5X, L5X, S5X;
};

struct CLK_FRAME {
	char mode[3];
	char name[5];
	int val_num;
	double CLK_OFF, CLK_DRI;
};

struct SP3_FRAME {
	double X, Y, Z, dt;
};

const double station_xyz[3] = { -1641945.2000, -3664804.1000, 4940009.3000 };
const double station_blh[3] = { 0.891513720275876, -1.99201147820819, 1118.80358502921};

const int num_of_date = 14;
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


UTC utc_obs = UTC::current_utc();
UTC utc_sp3 = UTC::current_utc();
UTC utc_clk = UTC::current_utc();

OBS_FRAME OBS[GPS_MAX_PRN];
SP3_FRAME SP3[GPS_MAX_PRN];
CLK_FRAME CLK[GPS_MAX_PRN];
CLK_FRAME REC;

bool obs_available[GPS_MAX_PRN] = { false };
bool sp3_available[GPS_MAX_PRN] = { false };
bool clk_available[GPS_MAX_PRN] = { false };
bool rec_available = false;

int sat_num = 0;

FILE * ofp, * sfp, * cfp;
char line_buffer[MAX_LENGTH_OF_LINE] = "";
char SVN[4] = "";
CLK_FRAME temp_clk;
double seconds = 0.0;
int prn_check = 0;

double klobuchar(double * sat_pos, int gps_sec, double * ea)
{
	return 0;
}


void fetch_clk(FILE * fp)
{
	fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
	sscanf(line_buffer, "%s%s%d%d%d%d%d%lf%d%lf%lf",
		&temp_clk.mode, &temp_clk.name,
		&utc_clk.year, &utc_clk.month, &utc_clk.date, &utc_clk.hour, &utc_clk.minute, &seconds,
		&temp_clk.val_num, &temp_clk.CLK_OFF, &temp_clk.CLK_DRI
	);
	utc_clk.year = utc_clk.year - 2000;
	utc_clk.sec = (int)round(seconds);


	int origin = ftell(fp);
	for (;; origin = ftell(fp)) {
		
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

	fseek(fp, origin, SEEK_SET);
}

void fetch_sp3(FILE * fp)
{
	fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
	sscanf(line_buffer, "*  %d %d %d %d %d %lf",
		&utc_sp3.year, &utc_sp3.month, &utc_sp3.date, &utc_sp3.hour, &utc_sp3.minute, &seconds);
	utc_sp3.year = utc_sp3.year - 2000;
	utc_sp3.sec = (int)round(seconds);

	while (!utc_sp3.equals(&utc_clk))
	{
		for (int i = 0; i < GPS_MAX_PRN; i++)
		{
			fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		}

		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		sscanf(line_buffer, "*  %d %d %d %d %d %lf",
			&utc_sp3.year, &utc_sp3.month, &utc_sp3.date, &utc_sp3.hour, &utc_sp3.minute, &seconds);
		utc_sp3.year = utc_sp3.year - 2000;
		utc_sp3.sec = (int)round(seconds);
	}

	for (int i = 0; i < GPS_MAX_PRN; i++)
	{
		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		double * arr = (double*)(SP3 + i);
		sscanf(line_buffer, "PG%d %lf %lf %lf %lf", 
			&prn_check, arr, arr + 1, arr + 2, arr + 3);
		sp3_available[i] = true;
	}
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

	while (!utc_obs.equals(&utc_sp3))
	{
		for (int i = 0; i < sat_num; i++)
		{
			fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		}

		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
		sscanf(line_buffer, "> %d %d %d %d %d %lf %d %d",
			&utc_obs.year, &utc_obs.month, &utc_obs.date, &utc_obs.hour, &utc_obs.minute, &seconds, &health, &sat_num);
		utc_obs.year = utc_obs.year - 2000;
		utc_obs.sec = (int)round(seconds);
	}

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

bool skip_sp3_header(FILE * fp)
{
	for (int i = 0; i < SP3_HEADER_LINE; i++)
	{
		if (feof(fp)) return false;
		fgets(line_buffer, MAX_LENGTH_OF_LINE, fp);
	}
	return true;
}

void overall_check()
{

}

void erace()
{
	memset(OBS, 0, sizeof(OBS_FRAME) * GPS_MAX_PRN);
	memset(SP3, 0, sizeof(SP3_FRAME) * GPS_MAX_PRN);
	memset(CLK, 0, sizeof(CLK_FRAME) * GPS_MAX_PRN);

	memset(obs_available, 0, sizeof(bool) * GPS_MAX_PRN);
	memset(sp3_available, 0, sizeof(bool) * GPS_MAX_PRN);
	memset(clk_available, 0, sizeof(bool) * GPS_MAX_PRN);
}

void execute(int index)
{
	ofp = fopen(obs_files[index], "r");
	sfp = fopen(sp3_files[index], "r");
	cfp = fopen(clk_files[index], "r");
	skip_obs_header(ofp);
	skip_sp3_header(sfp);
	skip_clk_header(cfp);

	while (!feof(ofp) && !feof(sfp) && !feof(cfp))
	{
		fetch_clk(cfp);
		fetch_sp3(sfp);
		fetch_obs(ofp);
		
		int i = 0;

	}

}

int main()
{
	for (int i = 0; i < num_of_date; i++)
	{
		erace();
		execute(i);
	}
}

