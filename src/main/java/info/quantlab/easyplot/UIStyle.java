package info.quantlab.easyplot;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.LegendTitle;
import java.awt.Color;
import java.awt.Font;

public class UIStyle {

    private final String fontName = "Arial";
    private final Font titleFont;
    private final Font axisLabelFont;
    private final Font annotationFont;
    private final Font legendFont;
    private final Font tickLabelFont;

    private final Color chartBackgroundPaint;
    private final Color plotBackgroundPaint;

    public UIStyle(Font titleFont, Font axisLabelFont, Font annotationFont, Font legendFont, Font tickLabelFont, Color chartBackgroundPaint, Color plotBackgroundPaint) {
        this.titleFont = titleFont;
        this.axisLabelFont = axisLabelFont;
        this.annotationFont = annotationFont;
        this.legendFont = legendFont;
        this.tickLabelFont = tickLabelFont;
        this.chartBackgroundPaint = chartBackgroundPaint;
        this.plotBackgroundPaint = plotBackgroundPaint;
    }

    public static UIStyle ScaledUIStyle(final double scale) {
        Font titleFont			= new Font(Font.SANS_SERIF, Font.PLAIN, (int)Math.round(10*scale));
        Font axisLabelFont		= new Font(Font.SANS_SERIF, Font.PLAIN, (int)Math.round(10*scale));
        Font annotationFont		= new Font(Font.SANS_SERIF, Font.PLAIN, (int)Math.round(8*scale));
        Font legendFont		    = new Font(Font.SANS_SERIF, Font.PLAIN, (int)Math.round(9*scale));
        Font tickLabelFont		= new Font(Font.SANS_SERIF, Font.PLAIN, (int)Math.round(9*scale));

        Color chartBackgroundPaint	= new java.awt.Color(250, 250, 250);
        //		chartBackgroundPaint	= new java.awt.Color(255, 255, 255);
        Color plotBackgroundPaint		= new java.awt.Color(255, 255, 255);
        return new UIStyle(titleFont, axisLabelFont, annotationFont, legendFont, tickLabelFont, chartBackgroundPaint, plotBackgroundPaint);
    }


    public static void applyStyleToChart(final JFreeChart chart) {
        new net.finmath.plots.jfreechart.StyleGuide(1).applyStyleToChart2(chart);
    }

    public static void applyStyleToXYPlot(final XYPlot xyPlot) {
        new net.finmath.plots.jfreechart.StyleGuide(1).applyStyleToXYPlot2(xyPlot);
    }

    public void applyStyleToChart2(final JFreeChart chart) {
        chart.setBackgroundPaint(chartBackgroundPaint);

        if(chart.getTitle() != null) {
            chart.getTitle().setFont(titleFont);
        }

        final LegendTitle legend = chart.getLegend();
        if(legend != null) {
            legend.setBackgroundPaint(chartBackgroundPaint);
            legend.setItemFont(legendFont);
        }

        final XYPlot xyPlot = chart.getXYPlot();
        if(xyPlot != null) {
            applyStyleToXYPlot2(xyPlot);
        }
    }

    public void applyStyleToXYPlot2(final XYPlot xyPlot) {
        for(int i=0; i < xyPlot.getDomainAxisCount(); i++) {
            final ValueAxis axis = xyPlot.getDomainAxis(i);
            if(axis != null) {
                axis.setTickLabelFont(tickLabelFont);
                axis.setLabelFont(axisLabelFont);
            }
        }
        xyPlot.setDomainAxisLocation(AxisLocation.BOTTOM_OR_LEFT);

        for(int i=0; i < xyPlot.getRangeAxisCount(); i++) {
            final ValueAxis axis = xyPlot.getRangeAxis(i);
            if(axis != null) {
                axis.setTickLabelFont(tickLabelFont);
                axis.setLabelFont(axisLabelFont);
            }
        }
        xyPlot.setRangeAxisLocation(AxisLocation.BOTTOM_OR_LEFT);
        xyPlot.setBackgroundPaint(plotBackgroundPaint);
    }

    /**
     * @return Font for axis labels.
     */
    public Font getAxisLabelFont() {
        return axisLabelFont;
    }

    /**
     * @return Font for axis tick labels.
     */
    public Font getTickLabelFont() {
        return tickLabelFont;
    }

    /**
     * @return Title font.
     */
    public Font getTitleFont() {
        return titleFont;
    }

}
