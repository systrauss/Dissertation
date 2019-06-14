#ifndef EVENTDATA_H
#define EVENTDATA_H

#include "TObject.h"
#include "TRandom3.h"

///Traditional data acquisiton setup at Univ. of Notre Dame 2013
/**This is a multipurpose setup used at ND for data acquistion. It consists
 * of a ADC in slot ## and a TDC in slot ##
 */ 
class eventData : public TObject
{
	public:
		
		// Timing Spectra
		UShort_t tdc[32]; //!temporary
		UShort_t sili_time[6];
		UShort_t pulser_time;
		UShort_t ge_time[2];
		UShort_t n1_time;
		UShort_t n2_time;
		UShort_t bgo6_time;
		UShort_t bgo1_time;
		UShort_t bgo5_time;
		UShort_t bgo4_time;
		UShort_t buncher_time;

		// Energy Spectra
		UShort_t adc[32]; //!temporary
		Float_t sili_en[6];
		UShort_t sili_en_raw[6];
		UShort_t ge_en_raw[2];
		Float_t ge_en[2];
		UShort_t n1_en[2];
		UShort_t n2_en[2];
		UShort_t pulser_en;
		UShort_t veto_timing_en;
		UShort_t bgo6_en;
		UShort_t bgo5_en;
		UShort_t bgo1_en;
		UShort_t bgo4_en;

	public:
		///Default Constructor.
		eventData();
		///Reset the current values.
		void Reset();
		///Record the value from the specified location.
		void SetValue(int crate, int slot, int ch, int value);
		///Calibrates the data from one event.
		void Calibrate();



	ClassDef(eventData,1);
};

#endif


