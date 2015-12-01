#ifndef CRC_32_HEADER_GUARD
#define CRC_32_HEADER_GUARD
#include <stdint.h> 
#include <cstddef>
uint32_t crc32(const void* data, std::size_t length, uint32_t previousCrc32 = 0);
#endif
