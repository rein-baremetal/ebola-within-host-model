package com.rptebov.app

/** Thin Kotlin wrapper over the C++ engine (libebov.so). */
object Engine {
    init { System.loadLibrary("ebov") }
    @JvmStatic external fun version(): String
    /** mode 0 = summary, 1 = trajectory table, 2 = SVG dashboard */
    @JvmStatic external fun run(params: String, mode: Int): String
    @JvmStatic external fun screen(params: String, cfg: String, prefix: String): String
}
