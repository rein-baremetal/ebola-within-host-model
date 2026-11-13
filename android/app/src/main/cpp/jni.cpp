// JNI bridge: text in, text out. Calls the engine's C API (see test_api.cpp).
#include <jni.h>
#include <string>
extern "C" {
const char* ebov_version();
char* ebov_run(const char*, int);
char* ebov_screen(const char*, const char*, const char*);
void ebov_free(char*);
}

static std::string str(JNIEnv* env, jstring s) {
    const char* c = env->GetStringUTFChars(s, nullptr);
    std::string r = c ? c : "";
    if (c) env->ReleaseStringUTFChars(s, c);
    return r;
}
static jstring take(JNIEnv* env, char* p) {
    jstring r = env->NewStringUTF(p ? p : "");
    ebov_free(p);
    return r;
}

extern "C" JNIEXPORT jstring JNICALL Java_com_rptebov_app_Engine_version(JNIEnv* env, jclass) {
    return env->NewStringUTF(ebov_version());
}
extern "C" JNIEXPORT jstring JNICALL Java_com_rptebov_app_Engine_run(JNIEnv* env, jclass, jstring params, jint mode) {
    return take(env, ebov_run(str(env, params).c_str(), mode));
}
extern "C" JNIEXPORT jstring JNICALL Java_com_rptebov_app_Engine_screen(JNIEnv* env, jclass, jstring params, jstring cfg, jstring prefix) {
    return take(env, ebov_screen(str(env, params).c_str(), str(env, cfg).c_str(), str(env, prefix).c_str()));
}
