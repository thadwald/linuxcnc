/********************************************************************
* Description: xyzac-switchkins.c
*   Derived from a work by Fred Proctor & Will Shackleford
*   and patch on forum:
*   https://forum.linuxcnc.org/10-advanced-configuration/31813-tcp-5-axis-kinematics?start=170#149538
*   based on work of
*   Rushabh Loladia: https://github.com/rushabhGH?tab=repositories (switchKins branch)
*
*  switchkins_type: 0 ==> trivkins, 1 ==> XYZAC
*
* License: GPL Version 2
* Copyright (c) 2009 All rights reserved.
*
********************************************************************/

#include "motion.h"
#include "hal.h"
#include "rtapi.h"     /* RTAPI realtime OS API */
#include "rtapi_app.h" /* RTAPI realtime module decls */
#include "rtapi_math.h"
#include "rtapi_string.h"
#include "kinematics.h"

// joint number assignments
static int JX = -1;
static int JY = -1;
static int JZ = -1;

static int JA = -1;
static int JC = -1;

struct haldata {
    hal_float_t *x_rot_point;
    hal_float_t *y_rot_point;
    hal_float_t *z_rot_point;
    hal_float_t *z_offset;
    hal_float_t *y_offset;
    hal_float_t *tool_offset;
    hal_bit_t   *kinstype_is_0;
    hal_bit_t   *kinstype_is_1;
} *haldata;

static hal_u32_t switchkins_type;

int kinematicsSwitchable() {return 1;}

int kinematicsSwitch(int new_switchkins_type)
{
    switchkins_type = new_switchkins_type;
    switch (switchkins_type) {
        case 0: rtapi_print_msg(RTAPI_MSG_INFO,
                "kinematicsSwitch:Trivkins\n");
                *haldata->kinstype_is_0 = 1;
                *haldata->kinstype_is_1 = 0;
                break;
        case 1: rtapi_print_msg(RTAPI_MSG_INFO,
                "kinematicsSwitch:XYZAC\n");
                *haldata->kinstype_is_1 = 1;
                *haldata->kinstype_is_0 = 0;
                break;
       default: rtapi_print_msg(RTAPI_MSG_ERR,
                "kinematicsSwitch:BAD VALUE <%d>\n",
                switchkins_type);
                *haldata->kinstype_is_1 = 0;
                *haldata->kinstype_is_0 = 0;
                return -1; // FAIL
    }

    return 0; // 0==> no error
}

static
int trivKinematicsForward(const double *joints,
                          EmcPose * pos,
                          const KINEMATICS_FORWARD_FLAGS * fflags,
                          KINEMATICS_INVERSE_FLAGS * iflags)
{
    pos->tran.x = joints[JX];
    pos->tran.y = joints[JY];
    pos->tran.z = joints[JZ];
    pos->a      = joints[JA];
    pos->c      = joints[JC];
    return 0;
}

static
int xyzacKinematicsForward(const double *joints,
                           EmcPose * pos,
                           const KINEMATICS_FORWARD_FLAGS * fflags,
                           KINEMATICS_INVERSE_FLAGS * iflags)
{

    double x_rot_point = *(haldata->x_rot_point);
    double y_rot_point = *(haldata->y_rot_point);
    double z_rot_point = *(haldata->z_rot_point);
    double          dt = *(haldata->tool_offset);
    double          dy = *(haldata->y_offset);
    double          dz = *(haldata->z_offset);
    double       a_rad = joints[JA]*TO_RAD;
    double       c_rad = joints[JC]*TO_RAD;

    dz = dz + dt;

    pos->tran.x = + cos(c_rad)              * (joints[JX]      - x_rot_point)
                  + sin(c_rad) * cos(a_rad) * (joints[JY] - dy - y_rot_point)
                  + sin(c_rad) * sin(a_rad) * (joints[JZ] - dz - z_rot_point )
                  + sin(c_rad) * dy
                  + x_rot_point;

    pos->tran.y = - sin(c_rad)              * (joints[JX]      - x_rot_point )
                  + cos(c_rad) * cos(a_rad) * (joints[JY] - dy - y_rot_point)
                  + cos(c_rad) * sin(a_rad) * (joints[JZ] - dz - z_rot_point)
                  + cos(c_rad) * dy
                  + y_rot_point;

    pos->tran.z = + 0
                  - sin(a_rad) * (joints[JY] - dy - y_rot_point)
                  + cos(a_rad) * (joints[JZ] - dz - z_rot_point)
                  + dz
                  + z_rot_point;

    pos->a = joints[JA];
    pos->c = joints[JC];

    pos->b = 0;
    pos->w = 0;
    pos->u = 0;
    pos->v = 0;

    return 0;
}


int kinematicsForward(const double *joints,
                      EmcPose * pos,
                      const KINEMATICS_FORWARD_FLAGS * fflags,
                      KINEMATICS_INVERSE_FLAGS * iflags)
{
    switch (switchkins_type) {
       case 1: return xyzacKinematicsForward(joints, pos, fflags, iflags);break;
      default: return  trivKinematicsForward(joints, pos, fflags, iflags);
    }

    return 0;
}

static
int trivKinematicsInverse(const EmcPose * pos,
                          double *joints,
                          const KINEMATICS_INVERSE_FLAGS * iflags,
                          KINEMATICS_FORWARD_FLAGS * fflags)
{
    joints[JX] = pos->tran.x;
    joints[JY] = pos->tran.y;
    joints[JZ] = pos->tran.z;
    joints[JA] = pos->a;
    joints[JC] = pos->c;
    return 0;
}

static
int xyzacKinematicsInverse(const EmcPose * pos,
                           double *joints,
                           const KINEMATICS_INVERSE_FLAGS * iflags,
                           KINEMATICS_FORWARD_FLAGS * fflags)
{
    double x_rot_point = *(haldata->x_rot_point);
    double y_rot_point = *(haldata->y_rot_point);
    double z_rot_point = *(haldata->z_rot_point);
    double         dy  = *(haldata->y_offset);
    double         dz  = *(haldata->z_offset);
    double         dt  = *(haldata->tool_offset);
    double      a_rad  = pos->a*TO_RAD;
    double      c_rad  = pos->c*TO_RAD;

    dz = dz + dt;

    joints[JX] = + cos(c_rad)              * (pos->tran.x - x_rot_point)
                 - sin(c_rad)              * (pos->tran.y - y_rot_point)
                 + x_rot_point;

    joints[JY] = + sin(c_rad) * cos(a_rad) * (pos->tran.x - x_rot_point)
                 + cos(c_rad) * cos(a_rad) * (pos->tran.y - y_rot_point)
                 -              sin(a_rad) * (pos->tran.z - z_rot_point)
                 -              cos(a_rad) * dy
                 +              sin(a_rad) * dz
                 + dy
                 + y_rot_point;

    joints[JZ] = + sin(c_rad) * sin(a_rad) * (pos->tran.x - x_rot_point)
                 + cos(c_rad) * sin(a_rad) * (pos->tran.y - y_rot_point)
                 +              cos(a_rad) * (pos->tran.z - z_rot_point)
                 -              sin(a_rad) * dy
                 -              cos(a_rad) * dz
                 + dz
                 + z_rot_point;


    joints[JA] = pos->a;
    joints[JC] = pos->c;

    return 0;
}

int kinematicsInverse(const EmcPose * pos,
                      double *joints,
                      const KINEMATICS_INVERSE_FLAGS * iflags,
                      KINEMATICS_FORWARD_FLAGS * fflags)
{
    switch (switchkins_type) {
       case 1: return xyzacKinematicsInverse(pos, joints, iflags, fflags);break;
      default: return  trivKinematicsInverse(pos, joints, iflags, fflags);
    }

    return 0;
}

int kinematicsHome(EmcPose * world,
                   double *joint,
                   KINEMATICS_FORWARD_FLAGS * fflags,
                   KINEMATICS_INVERSE_FLAGS * iflags)
{
    *fflags = 0;
    *iflags = 0;

    return kinematicsForward(joint, world, fflags, iflags);
}

KINEMATICS_TYPE kinematicsType()
{
    // both Forward and Inverse available for xyzac-trt and trivkins:
    return KINEMATICS_BOTH;
}

#define MAX_JOINTS 5
#define REQUIRED_COORDINATES "XYZAC"
static char *coordinates = REQUIRED_COORDINATES;
RTAPI_MP_STRING(coordinates, "Axes-to-joints-ordering");

EXPORT_SYMBOL(kinematicsSwitchable);
EXPORT_SYMBOL(kinematicsSwitch);
EXPORT_SYMBOL(kinematicsType);
EXPORT_SYMBOL(kinematicsForward);
EXPORT_SYMBOL(kinematicsInverse);
MODULE_LICENSE("GPL");

static int comp_id;
int rtapi_app_main(void)
{
static int axis_idx_for_jno[EMCMOT_MAX_JOINTS];
#define DISALLOW_DUPLICATES 0
    int res = 0;
    int jno;

    if (map_coordinates_to_jnumbers(coordinates,
                                    EMCMOT_MAX_JOINTS,
                                    DISALLOW_DUPLICATES,
                                    axis_idx_for_jno)) {
       return -1; //mapping failed
    }

    for (jno=0; jno<MAX_JOINTS; jno++) {
      if (axis_idx_for_jno[jno] == 0) {JX = jno;}
      if (axis_idx_for_jno[jno] == 1) {JY = jno;}
      if (axis_idx_for_jno[jno] == 2) {JZ = jno;}
      if (axis_idx_for_jno[jno] == 3) {JA = jno;}
      if (axis_idx_for_jno[jno] == 5) {JC = jno;}
    }
    if ( JX<0 || JY<0 || JZ<0 || JA<0 || JC<0 ) {
        rtapi_print_msg(RTAPI_MSG_ERR,
             "xyzac-trt-kins: required  coordinates:%s\n"
             "                specified coordinates:%s\n",
             REQUIRED_COORDINATES,coordinates);
        return -1;
    }
    rtapi_print("\nxyzac-trt-switchkins %s assigns:\n",coordinates);
    for (jno=0; jno<EMCMOT_MAX_JOINTS; jno++) {
        if (axis_idx_for_jno[jno] == -1) break; //fini
        rtapi_print("   Joint %d ==> Axis %c\n",
                   jno,*("XYZABCUVW"+axis_idx_for_jno[jno]));
    }

    comp_id = hal_init("xyzac-trt-switchkins");
    if(comp_id < 0) return comp_id;

    haldata = hal_malloc(sizeof(struct haldata));

    switchkins_type = 0; // startup with default type

    // conform to pin names for xyzac-trt-kins
    if((res = hal_pin_float_new("xyzac-trt-kins.x-rot-point",
              HAL_IN, &(haldata->x_rot_point), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("xyzac-trt-kins.y-rot-point",
              HAL_IN, &(haldata->y_rot_point), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("xyzac-trt-kins.z-rot-point",
              HAL_IN, &(haldata->z_rot_point), comp_id)) < 0) goto error;

    if((res = hal_pin_float_new("xyzac-trt-kins.y-offset",
              HAL_IN, &(haldata->y_offset), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("xyzac-trt-kins.z-offset",
              HAL_IN, &(haldata->z_offset), comp_id)) < 0) goto error;
    if((res = hal_pin_float_new("xyzac-trt-kins.tool-offset",
              HAL_IN, &(haldata->tool_offset), comp_id)) < 0) goto error;

    if((res = hal_pin_bit_new("kinstype.is-0",
              HAL_OUT, &(haldata->kinstype_is_0), comp_id)) < 0) goto error;
    if((res = hal_pin_bit_new("kinstype.is-1",
              HAL_OUT, &(haldata->kinstype_is_1), comp_id)) < 0) goto error;


    hal_ready(comp_id);
    return 0;

error:
    hal_exit(comp_id);
    return res;
}

void rtapi_app_exit(void) { hal_exit(comp_id); }
