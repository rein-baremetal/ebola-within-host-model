package com.rptebov.app

import android.app.Activity
import android.graphics.Typeface
import android.os.Bundle
import android.view.ViewGroup
import android.webkit.WebView
import android.widget.Button
import android.widget.EditText
import android.widget.LinearLayout
import android.widget.ScrollView
import android.widget.TextView
import java.io.File
import kotlin.concurrent.thread

/** Minimal front end: edit parameters, run the engine on the phone, see the summary and the SVG plot. */
class MainActivity : Activity() {
    private val defaults = """
        tip_dose = 1e9
        tip_time = 1
        c_T = 0.25
        alpha = 0.05
        beta_T = 2e-6
        V_ext = 2e-4
        ext_floor = 1
        t_end = 30
    """.trimIndent() + "\n"

    private lateinit var params: EditText
    private lateinit var out: TextView
    private lateinit var web: WebView

    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        val pad = (12 * resources.displayMetrics.density).toInt()
        val root = LinearLayout(this).apply { orientation = LinearLayout.VERTICAL; setPadding(pad, pad, pad, pad) }

        root.addView(TextView(this).apply { text = "RPT-EBOV"; textSize = 22f; setTypeface(null, Typeface.BOLD) })
        root.addView(TextView(this).apply {
            text = "Engine: ${Engine.version()}\nSimulated trajectories under stated assumptions; not medical advice."
            textSize = 12f
        })
        params = EditText(this).apply {
            setText(defaults); typeface = Typeface.MONOSPACE; textSize = 13f
            minLines = 8; gravity = android.view.Gravity.TOP
        }
        root.addView(params, ViewGroup.LayoutParams.MATCH_PARENT, ViewGroup.LayoutParams.WRAP_CONTENT)

        val row = LinearLayout(this).apply { orientation = LinearLayout.HORIZONTAL }
        fun button(label: String, action: () -> Unit) = Button(this).apply {
            text = label; setOnClickListener { action() }
        }.also { row.addView(it, LinearLayout.LayoutParams(0, ViewGroup.LayoutParams.WRAP_CONTENT, 1f)) }
        button("Run") { runEngine(0) }
        button("Plot") { runEngine(2) }
        button("Screen 200") { runScreen() }
        root.addView(row)

        out = TextView(this).apply { typeface = Typeface.MONOSPACE; textSize = 12f; setTextIsSelectable(true) }
        root.addView(out)
        web = WebView(this).apply {
            settings.useWideViewPort = true; settings.loadWithOverviewMode = true
            settings.builtInZoomControls = true; settings.displayZoomControls = false
        }
        root.addView(web, ViewGroup.LayoutParams.MATCH_PARENT, (420 * resources.displayMetrics.density).toInt())

        setContentView(ScrollView(this).apply { addView(root) })
    }

    private fun busy(msg: String) { out.text = msg }

    private fun runEngine(mode: Int) {
        val p = params.text.toString()
        busy("Running ...")
        thread {
            val res = Engine.run(p, mode)
            runOnUiThread {
                if (mode == 2 && res.startsWith("<")) {
                    out.text = "Plot ready (pinch to zoom)."
                    web.loadDataWithBaseURL(null, "<html><body style='margin:0'>$res</body></html>", "text/html", "utf-8", null)
                } else out.text = res
            }
        }
    }

    private fun runScreen() {
        val p = params.text.toString()
        busy("Screening 200 profiles ...")
        thread {
            val cfg = "samples = 200\nc_T = 0.1 5 log\ntip_dose = 1e6 1e9 log\ntip_time = 0 3\nbeta_T = 2.5e-7 2e-6 log\n"
            val prefix = File(filesDir, "screen").absolutePath
            val res = Engine.screen(p, cfg, prefix)
            runOnUiThread { out.text = res }
        }
    }
}
