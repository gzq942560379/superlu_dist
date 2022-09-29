#pragma once

#include <crts.h>
#include "sw/slave_param.h"

extern void SLAVE_FUN(remain_scatter_0_naive());
extern void SLAVE_FUN(remain_scatter_1_index_compress());
extern void SLAVE_FUN(remain_scatter_2_index_compress_dma());

extern void SLAVE_FUN(lookAhead_scatter_0_naive());
extern void SLAVE_FUN(lookAhead_scatter_1_index_compress());
extern void SLAVE_FUN(lookAhead_scatter_2_index_compress_dma());
