/*
******************************************************************************
File:     main.c
Info:     Generated by Atollic TrueSTUDIO(R) 7.1.1   2018-07-30

The MIT License (MIT)
Copyright (c) 2009-2017 Atollic AB

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

******************************************************************************
*/

/* Includes */
#include "stm32f4xx.h"
#include "tm_stm32f4_dac_signal.h"
#include "tm_stm32f4_dma.h"
#include "stm32f4xx_dac.h"
/* Private macro */
/* Private variables */
/* Private function prototypes */
/* Private functions */
uint16_t counter = 0;
uint8_t num_pulses = 16;
/**
**===========================================================================
**
**  Abstract: main program
**
**===========================================================================
*/
int main(void)
{
 /* Initialize system */
	SystemInit();

	/* Initialize DAC1, use TIM4 for signal generation */
	TM_DAC_SIGNAL_Init(TM_DAC1, TIM6);

	/* Output predefined sine signal with frequency of 40kHz */
	TM_DAC_SIGNAL_SetSignal(TM_DAC1, TM_DAC_SIGNAL_Signal_Sinus, 40000);

	while (1) {

		}

	}

