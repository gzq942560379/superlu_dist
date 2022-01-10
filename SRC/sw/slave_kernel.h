#pragma once

#include <athread.h>
#include "sw/slave_param.h"

extern void SLAVE_FUN(remain_scatter_0_naive());
extern void SLAVE_FUN(remain_scatter_1_ldm());
extern void SLAVE_FUN(remain_scatter_2_async());

extern void SLAVE_FUN(lookAhead_scatter_0_naive());
