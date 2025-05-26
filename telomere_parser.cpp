/* C++ module for regex search
   * This module provides a function to find the longest telomere repeat sequence
   * in a given DNA sequence, returning the type (C or G) and the sequence itself.
   * Author: Zhang Yuanfeng
   * Email: zhangyuanfeng1997@foxmail.com
   * Date: 2025-05-26
*/
#include <Python.h>
#include <string>
#include <regex>
#include <tuple>
#include <algorithm>

// 匹配结果结构体
struct TelomereMatch {
    std::string type;
    std::string sequence;
};


std::string find_longest_match(const std::string& seq, const std::regex& pattern) {
    std::sregex_iterator begin(seq.begin(), seq.end(), pattern), end;
    if (begin == end) return {};
    auto comp = [](const std::smatch& a, const std::smatch& b) {
        return a.length() < b.length();
    };

    return std::max_element(begin, end, comp)->str();
}


TelomereMatch telomere_strand_seq(const std::string& seq) {
    static const std::regex c_pattern(R"(([CT]{3}TAA)+)", std::regex::optimize);
    static const std::regex g_pattern(R"((TTA[GA]{3})+)", std::regex::optimize);

    // 查找最长匹配
    std::string longest_c = find_longest_match(seq, c_pattern);
    std::string longest_g = find_longest_match(seq, g_pattern);
    
    // 比较并返回结果
    if (!longest_c.empty() && (longest_c.length() >= longest_g.length() || longest_g.empty())) {
        return {"C", longest_c};
    } else if (!longest_g.empty()) {
        return {"G", longest_g};
    }
    return {"", ""};
}


// Python接口函数
static PyObject* py_telomere_strand_seq(PyObject* self, PyObject* args) {
    const char* seq_str;
    if (!PyArg_ParseTuple(args, "s", &seq_str)) {
        return NULL;
    }
    
    TelomereMatch result = telomere_strand_seq(seq_str);
    
    return Py_BuildValue("(ss)", result.type.c_str(), result.sequence.c_str());
}

// 方法表
static PyMethodDef StrandMethods[] = {
    {"telomere_strand_seq", py_telomere_strand_seq, METH_VARARGS, 
     "Find longest telomere repeat and type. Returns tuple (type, sequence)"},
    {NULL, NULL, 0, NULL}
};

// 模块定义
static struct PyModuleDef telomere_parser_module = {
    PyModuleDef_HEAD_INIT,
    "telomere_parser",
    "Telomere sequence parser module",
    -1,
    StrandMethods
};

// 初始化函数
PyMODINIT_FUNC PyInit_telomere_parser(void) {
    return PyModule_Create(&telomere_parser_module);
}
