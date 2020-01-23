#ifndef UNIT_TYPES_HEADER_GUARD
#define UNIT_TYPES_HEADER_GUARD
template <typename T, typename V> struct Unique
{
	V value;
	Unique(const V& value = V())
	:value(value)
	{}
	operator V() const
	{
		return value;
	}
};

//The marker pattern ID struct represents a unique identifier, assigned to a particular marker encoding.
struct markerPatternID_imp;
typedef Unique<markerPatternID_imp, int> markerPatternID;

//The funnel encoding struct represents a particular funnel, encoded as a bitfield
struct funnelEncoding_imp;
typedef Unique<funnelEncoding_imp, std::size_t> funnelEncoding;

//The funnel ID struct representse a unique identifier, assigned to a particular funnel encoding
struct funnelID_imp;
typedef Unique<funnelID_imp, int> funnelID;
#endif

