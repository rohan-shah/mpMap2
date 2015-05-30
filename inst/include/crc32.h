#ifndef CRC_32_HEADER_GUARD
#define CRC_32_HEADER_GUARD
#include <stdint.h> 
uint32_t crc32(const void* data, size_t length, uint32_t previousCrc32 = 0);
#endif