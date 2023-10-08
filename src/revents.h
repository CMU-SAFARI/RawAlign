#ifndef REVENTS_H
#define REVENTS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ri_events_s{
	uint32_t rid, length; //read id and length of the event values
	char *name; //name of the read

	float* values; //event  values of a read
} ri_events_t;

#ifdef __cplusplus
}
#endif
#endif //REVENTS_H
