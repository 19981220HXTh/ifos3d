/*
 * switch.h
 *
 *  Created on: Nov 21, 2017
 *      Author: bingo
 */

#ifndef SRC_SWITCH_H_
#define SRC_SWITCH_H_

#ifdef SW2
#define _MPE_
#else
#define _CPE_
#include <slave.h>
#endif




#endif /* SRC_SWITCH_H_ */
