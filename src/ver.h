#ifndef VER_H
#define VER_H

#ifdef __GNUC__
#define USED __attribute__((used))
#else
#define USED
#endif

#ifdef __cplusplus
extern "C" {
#endif

const char * const add_version_string(const char * const str);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#define VERSION(id) \
	static const char * const ver = add_version_string(id)
#else
#ifdef _MSC_VER
#define VERSION(id) \
#pragma comment(exestr, id)
#else
#define VERSION(id) \
	static const char * const USED ver = id
#endif
#endif

#endif /* VER_H */

