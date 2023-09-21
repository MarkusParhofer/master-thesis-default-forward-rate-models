/**
 * 
 */
package info.quantlab.easyplot;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Polygon;
import java.awt.Rectangle;
import java.awt.Shape;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jfree.chart.plot.DefaultDrawingSupplier;
import org.jfree.chart.plot.DrawingSupplier;

import net.finmath.plots.GraphStyle;
import net.finmath.plots.Named;
import net.finmath.plots.Plot2D;
import net.finmath.plots.Plotable2D;
import net.finmath.plots.PlotableFunction2D;
import net.finmath.plots.PlotablePoints2D;


/**
 * Does not work because of the set Method of List<...>.
 * @author Markus Parhofer
 *
 */
public class EasyPlot2D extends Plot2D {

	private static DrawingSupplier staticDrawingSupplier = new DefaultDrawingSupplier();
	
	private DrawingSupplier unstaticDrawingSupplier;
	
	private List<Plotable2D> plotables;
	
	// Multiple Plots
	public EasyPlot2D(final List<Plotable2D> plotables) {
		super(plotables);
		this.plotables = plotables;
		unstaticDrawingSupplier = staticDrawingSupplier;
		staticDrawingSupplier = new DefaultDrawingSupplier();
	}
	
	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final List<Named<DoubleUnaryOperator>> doubleUnaryOperators) {
		this(doubleUnaryOperators.stream().map(namedFunction -> 
		{ return new PlotableFunction2D(xmin, xmax, numberOfPointsX, namedFunction, 
				getDefaultGraphStyle(doubleUnaryOperators.indexOf(namedFunction), true, staticDrawingSupplier)); })
				.collect(Collectors.toList()));
	}

	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final DoubleUnaryOperator[] doubleUnaryOperators) {
		this(xmin, xmax, numberOfPointsX,
				Arrays.stream(doubleUnaryOperators).map(operator -> new Named<DoubleUnaryOperator>("", operator)).collect(Collectors.toList())
				);
	}

	
	// Single Plots
	public EasyPlot2D(final Plotable2D plotable) {
		this(Collections.singletonList(plotable));
	}
	
	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final Named<DoubleUnaryOperator> function, final GraphStyle style) {
		this(new PlotableFunction2D(xmin, xmax, numberOfPointsX, function, style));
	}
	
	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final Named<DoubleUnaryOperator> function) {
		this(new PlotableFunction2D(xmin, xmax, numberOfPointsX, function, getDefaultGraphStyle(0, true, staticDrawingSupplier)));
	}
	
	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final DoubleUnaryOperator function, final GraphStyle style) {
		this(new PlotableFunction2D(xmin, xmax, numberOfPointsX, new Named<DoubleUnaryOperator>("", function), style));
	}

	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final DoubleUnaryOperator function) {
		this(xmin, xmax, numberOfPointsX, Collections.singletonList(new Named<DoubleUnaryOperator>("",function)));
	}
	
	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final Function<Double, Double> function) {
		this(xmin, xmax, numberOfPointsX, Collections.singletonList(new Named<DoubleUnaryOperator>("", new DoubleUnaryOperator() {
			@Override
			public double applyAsDouble(double operand) {
				return function.apply(operand);
			}
		})));
	}
	
	
	// Plotting Points
	public EasyPlot2D(String name, final double[] xPoints, final double[] yPoints, final GraphStyle style) {
		this(PlotablePoints2D.of(name, xPoints, yPoints, style));
	}
	
	public EasyPlot2D(String name, final double[] xPoints, final double[] yPoints) {
		this(PlotablePoints2D.of(name, xPoints, yPoints, getDefaultGraphStyle(0, false, staticDrawingSupplier)));
	}
	
	public EasyPlot2D(final double[] xPoints, final double[] yPoints) {
		this(PlotablePoints2D.of("", xPoints, yPoints, getDefaultGraphStyle(0, false, staticDrawingSupplier)));
	}

	public EasyPlot2D addPlot(final Plotable2D plotable) {
		Stream<Plotable2D> stringStream = Stream.concat(plotables.stream(), Stream.of(plotable));
		plotables = stringStream.collect(Collectors.toList());
		update();
		return this;
	}
	
	public EasyPlot2D addPlot(final double xmin, final double xmax, final int numberOfPointsX, final Named<DoubleUnaryOperator> function, final GraphStyle style) {
		return addPlot(new PlotableFunction2D(xmin, xmax, numberOfPointsX, function, style));
	}
	
	public EasyPlot2D addPlot(final double xmin, final double xmax, final int numberOfPointsX, final Named<DoubleUnaryOperator> function) {
		return addPlot(new PlotableFunction2D(xmin, xmax, numberOfPointsX, function, getDefaultGraphStyle(plotables.size(), true, unstaticDrawingSupplier)));
	}
	
	@Override
	public Plot2D update(final List<Plotable2D> plotables) {
		this.plotables = plotables;
		return super.update(plotables);
	}
	
	public void update() {
		super.update(plotables);
	}
	
	public GraphStyle changePlotStyle(final int plotIndex, final GraphStyle newStyle) {
		Plotable2D plotToChange = plotables.get(plotIndex);
		
		final boolean styleWasSet = plotToChange.getStyle() != null;
		GraphStyle returnStyle = styleWasSet ? plotToChange.getStyle() : getDefaultGraphStyle(plotIndex, true, unstaticDrawingSupplier);
		
		final double xMin = plotToChange.getSeries().get(0).getX();
		final double xMax = plotToChange.getSeries().get(plotToChange.getSeries().size() - 1).getX();
		
		Named<DoubleUnaryOperator> newFunction = new Named<DoubleUnaryOperator>(plotToChange.getName(), 
				new DoubleUnaryOperator() {
			@Override
			public double applyAsDouble(double operand) {
				for(int index = 0; index < plotToChange.getSeries().size(); index++) {
					if(plotToChange.getSeries().get(index).getX() == operand)
						return plotToChange.getSeries().get(index).getY();
				}
				return 0;
			}
		});
		
		Plotable2D newPlot = new PlotableFunction2D(xMin, xMax, plotToChange.getSeries().size(), newFunction, newStyle);
		List<Plotable2D> newPlotables = plotables.stream().map(plotable -> {
			if(plotable == plotToChange)
				return newPlot;
			else
				return plotable;
		}).collect(Collectors.toList());
		plotables = newPlotables;
		update();
		return returnStyle;
	}
	
	public Color changePlotColor(final int plotIndex, final Color newColor) {
		
		Plotable2D plotToChange = plotables.get(plotIndex);
		
		final boolean styleWasSet = plotToChange.getStyle() != null;
		
		GraphStyle newGraphStyle = new GraphStyle(styleWasSet? plotToChange.getStyle().getShape() : null, styleWasSet? plotToChange.getStyle().getStroke() : new BasicStroke(2.0f), newColor);
		
		GraphStyle oldStyle = changePlotStyle(plotIndex, newGraphStyle);
		return oldStyle.getColor();
	}
	
	public Color changePlotColor(final Color newColor) {
		return changePlotColor(0, newColor);
	}
	
	public static Color getDefaultColor(final int colorIndex) {
		switch (colorIndex % 10) {
		case 0:
			return new java.awt.Color(0, (int)(0.4470*255), (int)(0.7410*255));
		case 1:
			return new java.awt.Color((int)(0.8500*255), (int)(0.3250*255), (int)(0.0980*255));
		case 2:
			return new java.awt.Color((int)(0.9290*255), (int)(0.6940*255), (int)(0.1250*255));
		case 3:
			return new java.awt.Color((int)(0.4940*255), (int)(0.1840*255), (int)(0.5560*255));
		case 4:
			return new java.awt.Color((int)(0.4660*255), (int)(0.6740*255), (int)(0.1880*255));
		case 5:
			return new java.awt.Color((int)(0.3010*255), (int)(0.7450*255), (int)(0.9330*255));
		case 6:
			return new java.awt.Color((int)(0.6350*255), (int)(0.0780*255), (int)(0.1840*255));
		case 7:
			return new java.awt.Color(255, 0, 0);
		case 8:
			return new java.awt.Color(0, 255, 0);
		case 9:
			return new java.awt.Color(0, 0, 255);
		default:
			return new java.awt.Color(0, 0,  0);
		}
		
	}
	
	public static GraphStyle getDefaultGraphStyle(final int funcIndex, final boolean asLine, DrawingSupplier drawingSupplier) {
		if(asLine) {
			return new GraphStyle(drawingSupplier.getNextShape(), new BasicStroke(2.0f), getDefaultColor(funcIndex));
		} else {
			return new GraphStyle(drawingSupplier.getNextShape(), null, getDefaultColor(funcIndex * 10 + (((int)Math.floor(funcIndex / 10.0)) % 10)));
		}
	}
	
	public static Color getColor(final double red, final double green, final double blue) {
		return new java.awt.Color((int)(red*255), (int)(green*255), (int)(blue*255));
	}
	
	public static Color getColor(final int red, final int green, final int blue) {
		return new java.awt.Color(red, green, blue);
	}
	
}
