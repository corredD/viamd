#include "string_utils.h"
#include <core/common.h>
#include <fstream>
#include <ctype.h>

#ifdef WIN32
#pragma warning(disable:4996) // Unsafe strncpy
#endif

#define MIN(x,y) ((x < y) ? (x) : (y))
#define MAX(x,y) ((x > y) ? (x) : (y))

static bool internal_compare(const char* str_a, const char* str_b, int64 len, bool ignore_case) {
	if (ignore_case) {
		for (int64 i = 0; i < len; i++) {
			if (tolower(str_a[i]) != tolower(str_b[i])) return false;
		}
	}
	else {
		for (int64 i = 0; i < len; i++) {
			if (str_a[i] != str_b[i]) return false;
		}
	}
	return true;
}

bool compare(CString str_a, CString str_b, bool ignore_case) {
	int64 len = MIN(str_a.count, str_b.count);
	return internal_compare(str_a, str_b, len, ignore_case);
}

bool compare_n(CString str_a, CString str_b, int64 n, bool ignore_case) {
	int64 len = MIN(str_a.count, MIN(str_b.count, n));
	return internal_compare(str_a, str_b, len, ignore_case);
}

void copy(String dst, CString src) {
	ASSERT(dst.data != 0);
	ASSERT(src.data != 0);
	ASSERT(dst.count >= src.count);
	// Should one use strncpy or memcpy?
	//strncpy(dst.data, src.data, src.count);
	memcpy(dst.data, src.data, src.count);
}

bool get_line(String& line, CString& str) {
	const char* str_beg = str.data;
	const char* str_end = str.data + str.count;

	if (str_beg == str_end) {
		line = {};
		return false;
	}

	const char* line_beg = str_beg;
	const char* line_end = line_beg;

	while (line_end < str_end && *line_end != '\n') ++line_end;

	str_beg = MIN(line_end + 1, str_end);

	// If we have a preceding '\r' we need to remove that as well
	if (line_end > line_beg && *(line_end - 1) == '\r') {
		--line_end;
	}

	line.count = line_end - line_beg;
	line.data = (char*)memcpy(line.data, line_beg, line.count);

	str.data = str_beg;
	str.count = str_end - str_beg;

	return true;
}

float to_float(CString str) {
	// Make sure that the string passed into atof is zero-terminated
	char buffer[64] = {};
	memcpy(buffer, str.data, MIN(63, str.count));
	return (float)atof(buffer);
}

int to_int(CString str) {
	// Make sure that the string passed into atoi is zero-terminated
	char buffer[64] = {};
	memcpy(buffer, str.data, MIN(63, str.count));
	return atoi(buffer);
}

CString trim(CString str) {
	const char* beg = str.data;
	const char* end = str.data + str.count;

	while (beg < end && isspace(*beg)) ++beg;
	while (end > beg && isspace(*(end-1))) --end;

	return CString(beg, end - beg);
}

String trim(String str) {
	char* beg = str.data;
	char* end = str.data + str.count;

	while (beg < end && isspace(*beg)) ++beg;
	while (end > beg && isspace(*(end-1))) --end;

	return String(beg, end - beg);
}

String read_textfile(CString filename, Allocator& alloc) {
	std::ifstream file(filename);
	if (!file) return {};

	file.seekg(0, std::ios::end);
	int64 size = file.tellg();
	file.seekg(0, std::ios::beg);

	char* data = (char*)alloc.alloc(size);
	file.read(data, size);

	return {data, size};
}

CString get_directory(CString url) {
    if (url.count == 0) {
        return url;
    }
    
    url = trim(url);
    
    const char* beg = url.begin();
    const char* end = url.end() - 1;
    
    while (end != beg && *end != '\\' && *end != '/') {
        end--;
    }
    
    return CString(beg, end - beg);
}

CString get_file(CString url) {
    if (url.count == 0) {
        return url;
    }

    url = trim(url);
    
    const char* beg = url.end() - 1;
    const char* end = url.end();
    
    while (beg != url.begin() && *beg != '\\' && *beg != '/') {
        beg--;
    }
    
    return CString(beg, end - beg);
}

CString get_file_without_extension(CString url) {
    if (url.count == 0) {
        return url;
    }
    
    url = trim(url);
    
    const char* beg = url.end() - 1;
    const char* end = url.end();
    
    while (beg != url.begin() && *beg != '\\' && *beg != '/') beg--;
    while (end != beg && *beg != '.') end--;
    
    return CString(beg, end - beg);
}

CString get_file_extension(CString url) {
    if (url.count == 0) {
        return url;
    }
    
    url = trim(url);
    
    const char* beg = url.end() - 1;
    const char* end = url.end();
    
    while (beg != url.begin() && *beg != '.') beg--;
    
    if (beg == url.begin()) {
        return CString();
    }
    
    return CString(beg, end - beg);
}