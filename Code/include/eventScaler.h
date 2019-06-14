#ifndef EVENTSCALER_H
#define EVENTSCALER_H

#define NUM_OF_SCALER_CH 26

#include "TObject.h"

class eventScaler : public TObject
{
	public:
		UShort_t scalerStartTime;
		UShort_t scalerEndTime;
		UShort_t scalers[NUM_OF_SCALER_CH];

		eventScaler();
		void Reset();
		void SetEndTime(unsigned int time);
		void SetStartTime(unsigned int time);
		void SetValue(int channel, unsigned int value);

	ClassDef(eventScaler,1)
};

#endif
