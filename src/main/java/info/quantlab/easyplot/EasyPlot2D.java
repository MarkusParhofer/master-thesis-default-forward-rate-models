/**
 * 
 */
package info.quantlab.easyplot;
import java.awt.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import info.quantlab.debug.Debug;
import net.finmath.plots.*;
import net.finmath.plots.axis.NumberAxis;
import org.jfree.chart.plot.DefaultDrawingSupplier;
import org.jfree.chart.plot.DrawingSupplier;


/**
 * @author Markus Parhofer
 *
 */
public class EasyPlot2D extends Plot2D {

	private static DrawingSupplier staticDrawingSupplier = new DefaultDrawingSupplier();

	private final DrawingSupplier unstaticDrawingSupplier;

	private static Shape lastShape;

	private String title = null;
	private String xAxisName = null;
	private String yAxisName = null;
	private boolean legend = false;
	private double xMin = Double.NaN, xMax = Double.NaN;
	private double yMin = Double.NaN, yMax = Double.NaN;

	private List<Plotable2D> plotables;

	// Multiple Plots
	public EasyPlot2D(final List<Plotable2D> plotables) {
		super(plotables);
		if (plotables.get(0) instanceof Plotable2DCloneable)
			this.plotables = plotables;
		else {
			update(plotables);
		}
		unstaticDrawingSupplier = staticDrawingSupplier;
		staticDrawingSupplier = new DefaultDrawingSupplier();

	}

	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final List<Named<DoubleUnaryOperator>> doubleUnaryOperators) {
		this(doubleUnaryOperators.stream().map(namedFunction ->
						new PlotableFunction2D(xmin, xmax, numberOfPointsX, namedFunction,
								getDefaultGraphStyle(doubleUnaryOperators.indexOf(namedFunction), true, staticDrawingSupplier)))
				.collect(Collectors.toList()));
	}

	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final DoubleUnaryOperator[] doubleUnaryOperators) {
		this(xmin, xmax, numberOfPointsX,
				Arrays.stream(doubleUnaryOperators).map(operator -> new Named<>("", operator)).collect(Collectors.toList())
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
		this(new PlotableFunction2D(xmin, xmax, numberOfPointsX, new Named<>("", function), style));
	}

	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final DoubleUnaryOperator function) {
		this(xmin, xmax, numberOfPointsX, Collections.singletonList(new Named<>("", function)));
	}

	public EasyPlot2D(final double xmin, final double xmax, final int numberOfPointsX, final Function<Double, Double> function) {
		this(xmin, xmax, numberOfPointsX, Collections.singletonList(new Named<>("", function::apply)));
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

	public EasyPlot2D(final double[] xPoints, final double[][] yPoints) {
		this(IntStream.range(0, yPoints.length).mapToObj(i ->
				PlotablePoints2D.of("Series " + i, xPoints, yPoints[i],
						getDefaultGraphStyle(i, false, staticDrawingSupplier))).collect(Collectors.toList()));
	}

	public EasyPlot2D(final String[] plotNames, final double[] xPoints, final double[][] yPoints) {
		this(IntStream.range(0, yPoints.length).mapToObj(i ->
				PlotablePoints2D.of(plotNames[i], xPoints, yPoints[i],
						getDefaultGraphStyle(i, false, staticDrawingSupplier))).collect(Collectors.toList()));
	}

	public EasyPlot2D addPlot(final Plotable2D plotable) {
		Stream<Plotable2D> plotableStream = Stream.concat(plotables.stream(), Stream.of(new Plotable2DCloneable(plotable)));
		plotables = plotableStream.collect(Collectors.toList());
		update();
		return this;
	}

	public EasyPlot2D addPlot(final List<Plotable2D> newPlotables) {
		List<Plotable2D> plotable2DS = newPlotables.get(0) instanceof Plotable2DCloneable ? newPlotables :
				Plotable2DCloneable.listOf(newPlotables, unstaticDrawingSupplier);
		Stream<Plotable2D> plotStream = Stream.concat(plotables.stream(), plotable2DS.stream());
		plotables = plotStream.collect(Collectors.toList());
		update();
		return this;
	}

	public EasyPlot2D addPlot(final double xmin, final double xmax, final int numberOfPointsX, final Named<DoubleUnaryOperator> function, final GraphStyle style) {
		return addPlot(new PlotableFunction2D(xmin, xmax, numberOfPointsX, function, style));
	}

	public EasyPlot2D addPlot(final double xmin, final double xmax, final int numberOfPointsX, final Named<DoubleUnaryOperator> function) {
		return addPlot(new PlotableFunction2D(xmin, xmax, numberOfPointsX, function, getDefaultGraphStyle(plotables.size(), true, unstaticDrawingSupplier)));
	}

	public EasyPlot2D addPlot(final double xmin, final double xmax, final int numberOfPointsX, List<Named<DoubleUnaryOperator>> functions) {
		List<Plotable2D> newPlots = functions.stream().map(operator -> (Plotable2D) (new PlotableFunction2D(
				xmin, xmax, numberOfPointsX, operator, getDefaultGraphStyle(functions.indexOf(operator), true, unstaticDrawingSupplier)))).toList();
		return addPlot(newPlots);
	}

	public EasyPlot2D addPlot(String name, final double[] xPoints, final double[] yPoints, final GraphStyle style) {
		return addPlot(PlotablePoints2D.of(name, xPoints, yPoints, style));
	}

	public EasyPlot2D addPlot(String name, final double[] xPoints, final double[] yPoints) {
		return addPlot(name, xPoints, yPoints, getDefaultGraphStyle(getNumberOfPlots(), false, unstaticDrawingSupplier));
	}

	public EasyPlot2D addPlot(final double[] xPoints, final double[] yPoints) {
		return addPlot("", xPoints, yPoints);
	}

	public EasyPlot2D removePlot(final int index) {
		if (index < 0 || index >= getNumberOfPlots()) {
			return this;
		}
		ArrayList<Plotable2D> plots = new ArrayList<>(plotables);
		plots.remove(index);
		plotables = plots;
		update();
		return this;
	}

	public boolean writeRecreationJavaFile() {
		String timeStamp = new SimpleDateFormat("MM-dd-yyyy").format(new Date());
		timeStamp = timeStamp.replaceAll("[^abcdefghjklmnopqrstuvwxyzABCDEFGHJKLMNOPQRSTUVWXYZ1234567890_]", "_");
		String pathString = "src\\main\\java\\info\\quantlab\\easyplot\\recreatefiles\\";
		File directory = new File(pathString);
		boolean directoryExists = directory.exists();
		if(!directoryExists)
			directoryExists = directory.mkdirs();

		if(directoryExists) {
			String className = title.replaceAll("[^abcdefghjklmnopqrstuvwxyzABCDEFGHJKLMNOPQRSTUVWXYZ1234567890_]", "_");
			className += timeStamp;
			File file = new File(pathString + className + ".java");
		boolean result;
		try {
			result = file.createNewFile();
		} catch (IOException ioEx) {
			Debug.logln(ioEx.getMessage());
			return false;
		} catch (SecurityException secEx) {
			Debug.logln(secEx.getMessage());
			Debug.logln("Access to File Denied!");
			return false;
		}
		if (!result)
			return false;

		try {
			FileWriter myWriter = new FileWriter(file);
			// Package
			{
				String packageName = "package info.quantlab.easyplot.recreatefiles;\n";
				myWriter.write(packageName);
			}
			// Imports
			{
				String imp = "import info.quantlab.easyplot.EasyPlot2D;\n" +
						"import net.finmath.plots.GraphStyle;\n" +
						"import net.finmath.plots.Plotable2D;\n" +
						"import net.finmath.plots.PlotablePoints2D;\n" +
						"import org.jfree.chart.plot.DefaultDrawingSupplier;\n" +
						"import org.jfree.chart.plot.DrawingSupplier;\n" +
						"import java.util.List;\n" +
						"import java.util.ArrayList;\n\n\n";
				myWriter.write(imp);
			}
			// Class name and main function and List of plotables
			{
				String classN = "public class " + className + " {\n\n" +
						"public static void main(String[] args) {\n\n" +
						"DrawingSupplier drawingSupplier = new DefaultDrawingSupplier();\n" +
						"ArrayList<Plotable2D> plots = new ArrayList<>();\n\n";
				myWriter.write(classN);
			}
			// Write Plotables
			for(int i=0; i < getNumberOfPlots(); i++) {
				Plotable2DCloneable p;
				try{
					p = (Plotable2DCloneable)plotables.get(i);
				} catch (ClassCastException ex) {
					System.out.println("Plot is not of Type \"Plotable2DCloneable\". There must be a bug somewhere.");
					p = Plotable2DCloneable.of(plotables.get(i));
					update(plotables);
				}
				myWriter.write("{\n");
				// Points:
				{
					final List<Point2D> series = p.getSeries();
					StringBuilder pointsX = new StringBuilder("final double[] xPoints = new double[] {\n");
					StringBuilder pointsY = new StringBuilder("final double[] yPoints = new double[] {\n");
					for (int j = 0; j < series.size(); j++) {
						String div = j == series.size() - 1 ? "};" : ",\n";
						pointsX.append(String.format(Locale.ENGLISH, "%30.15f%2s", series.get(j).getX(), div));
						pointsY.append(String.format(Locale.ENGLISH, "%30.15f%2s", series.get(j).getY(), div));
					}
					myWriter.write(pointsX.toString() + "\n");
					myWriter.write(pointsY.toString() + "\n");
				}

				// Style
				{
					final GraphStyle style = p.getStyle();
					String asLine = p.getKind() == Plotable2DCloneable.Kind.FUNCTION? "true" : "false";
					String sStyle = "final GraphStyle style = EasyPlot2D.getDefaultGraphStyle(" + i +", "+ asLine + ", drawingSupplier);\n";
					String comment = "/* Original Style: " + style.toString() + "*/\n";
					myWriter.write(sStyle);
					myWriter.write(comment);
				}
				String plot = "final Plotable2D p = PlotablePoints2D.of(\"" +
						p.getName() + "\", xPoints, yPoints, style);\n";
				myWriter.write(plot);
				String attachToList = "plots.add(p);\n";
				myWriter.write(attachToList);
				myWriter.write("}\n");
			}
			myWriter.write("EasyPlot2D myPlot = new EasyPlot2D(plots);\n");
			if(title != null) {
				myWriter.write("myPlot.setTitle(\"" + title + "\");\n");
			}
			if(xAxisName != null) {
				myWriter.write("myPlot.setXAxisLabel(\"" + xAxisName + "\");\n");
			}
			if(yAxisName != null) {
				myWriter.write("myPlot.setYAxisLabel(\"" + yAxisName + "\");\n");
			}
			if(!Double.isNaN(xMin) && !Double.isNaN(xMax)) {
				myWriter.write("myPlot.setXRange(" + xMin + ", " + xMax + ");\n");
			}
			if(!Double.isNaN(yMin) && !Double.isNaN(yMax)) {
				myWriter.write("myPlot.setYRange(" + yMin + ", " + yMax + ");\n");
			}
			myWriter.write("myPlot.setIsLegendVisible(" + (legend? "true" : "false") + ");\n");
			myWriter.write("myPlot.show();\n");
			String comment = "/* Original Plot: \"" + super.toString() + "\"*/\n";
			myWriter.write(comment);
			// close main
			myWriter.write("}\n");
			// close class
			myWriter.write("}\n");

			myWriter.close();
			System.out.println("Successfully wrote to the file.");
		} catch(IOException e)
		{
			System.out.println("An error occurred.");
			return false;
		}
		return true;
		}

		System.out.println("Cannot create directory for Plot");
		return false;
	}

	public void printPoints(int index) {
		final List<Point2D> points = plotables.get(index).getSeries();
		System.out.print("\nPlotable " + index + " x: {");
		for(int i=0; i < points.size(); i++) {
			String divider = i < points.size()-1 ? "," : "}";
			System.out.printf(" %20.11f%1s", points.get(i).getX(), divider);
		}
		System.out.print("\nPlotable " + index + " y: {");
		for(int i=0; i < points.size(); i++) {
			String divider = i < points.size()-1 ? "," : "}";
			System.out.printf(" %20.11f%1s", points.get(i).getY(), divider);
		}
		System.out.println();
	}

	@Override
	public Plot2D update(final List<Plotable2D> plotables) {
		this.plotables = Plotable2DCloneable.listOf(plotables, unstaticDrawingSupplier);
		return super.update(this.plotables);
	}

	public Plot2D update() {
		return super.update(new ArrayList<>(plotables));
	}

	public String changePlotName(final int plotIndex, final String name) {
		Plotable2DCloneable plotToChange = null;
		boolean error = false;
		try {
			plotToChange = (Plotable2DCloneable)plotables.get(plotIndex);
		} catch (ClassCastException ex) {
			System.out.println("Plot is not of Type \"Plotable2DCloneable\". There must be abug somewhere.");
			plotToChange = Plotable2DCloneable.of(plotables.get(plotIndex));
			error = true;
		}
		String returnName = plotToChange.getName();
		final Plotable2DCloneable plotToChangeC = plotToChange.getCloneWithModifiedName(name);

		plotables = plotables.stream().map(plotable -> {
			if(plotable == plotables.get(plotIndex))
				return plotToChangeC;
			else
				return plotable;
		}).collect(Collectors.toList());
		if(error) {
			update(plotables);
		} else {
			update();
		}
		return returnName;
	}

	public String changePlotName(final String name) {
		return changePlotName(0, name);
	}

	public GraphStyle changePlotStyle(final int plotIndex, final GraphStyle newStyle) {
		Plotable2DCloneable plotToChange = null;
		boolean error = false;
		try {
			plotToChange = (Plotable2DCloneable)plotables.get(plotIndex);
		} catch (ClassCastException ex) {
			System.out.println("Plot is not of Type \"Plotable2DCloneable\". There must be abug somewhere.");
			plotToChange = Plotable2DCloneable.of(plotables.get(plotIndex));
			error = true;
		}
		final boolean styleWasSet = plotToChange.getStyle() != null;
		GraphStyle returnStyle = styleWasSet ? plotToChange.getStyle() :
				getDefaultGraphStyle(plotIndex,
                        plotToChange.getKind() == Plotable2DCloneable.Kind.FUNCTION, unstaticDrawingSupplier);
		final Plotable2DCloneable plotToChangeC = plotToChange.getCloneWithModifiedStyle(newStyle);

		plotables = plotables.stream().map(plotable -> {
			if(plotable == plotables.get(plotIndex))
				return plotToChangeC;
			else
				return plotable;
		}).collect(Collectors.toList());
		if(error) {
			update(plotables);
		} else {
			update();
		}
		return returnStyle;
	}

	public int getNumberOfPlots() {
		return plotables.size();
	}

	public Color changePlotColor(final int plotIndex, Color newColor) {

		if(newColor == null) {
			newColor = getDefaultColor(plotIndex);
		}
		Plotable2D plotToChange = plotables.get(plotIndex);

		GraphStyle newGraphStyle = new GraphStyle(plotToChange.getStyle().getShape(), plotToChange.getStyle().getStroke(), newColor);
		
		GraphStyle oldStyle = changePlotStyle(plotIndex, newGraphStyle);
		return oldStyle.getColor();
	}
	
	public Color changePlotColor(final Color newColor) {
		return changePlotColor(0, newColor);
	}

	public Shape changePlotMarker(final int plotIndex, final Shape newShape) {

		Plotable2D plotToChange = plotables.get(plotIndex);
		Stroke stroke = plotToChange.getStyle().getStroke();
		if(stroke == null && newShape == null) {
			stroke = new BasicStroke(2.0f);
		}

		GraphStyle newGraphStyle = new GraphStyle(newShape, stroke, plotToChange.getStyle().getColor());

		GraphStyle oldStyle = changePlotStyle(plotIndex, newGraphStyle);

		return oldStyle.getShape();
	}

	public Stroke changePlotStroke(final int plotIndex, final Stroke newStroke) {

		Plotable2D plotToChange = plotables.get(plotIndex);
		Shape shape = plotToChange.getStyle().getShape();
		if(newStroke == null && shape == null) {
			shape = skipShape();
		}


		GraphStyle newGraphStyle = new GraphStyle(shape, newStroke, plotToChange.getStyle().getColor());

		GraphStyle oldStyle = changePlotStyle(plotIndex, newGraphStyle);

		return oldStyle.getStroke();
	}

	public Stroke changePlotStroke(final Stroke newStroke) {
		return changePlotStroke(0, newStroke);
	}

	public Shape changePlotMarker(final Shape newShape) {
		return changePlotMarker(0, newShape);
	}

	public Shape skipShape() {
		return unstaticDrawingSupplier.getNextShape();
	}

	public static Color getDefaultColor(final int colorIndex) {
        return switch (colorIndex % 10) {
            case 0 -> new Color(0, (int) (0.4470 * 255), (int) (0.7410 * 255));
            case 1 -> new Color((int) (0.8500 * 255), (int) (0.3250 * 255), (int) (0.0980 * 255));
            case 2 -> new Color((int) (0.9290 * 255), (int) (0.6940 * 255), (int) (0.1250 * 255));
            case 3 -> new Color((int) (0.4940 * 255), (int) (0.1840 * 255), (int) (0.5560 * 255));
            case 4 -> new Color((int) (0.4660 * 255), (int) (0.6740 * 255), (int) (0.1880 * 255));
            case 5 -> new Color((int) (0.3010 * 255), (int) (0.7450 * 255), (int) (0.9330 * 255));
            case 6 -> new Color((int) (0.6350 * 255), (int) (0.0780 * 255), (int) (0.1840 * 255));
            case 7 -> new Color(255, 0, 0);
            case 8 -> new Color(0, 255, 0);
            case 9 -> new Color(0, 0, 255);
            default -> new Color(0, 0, 0);
        };
		
	}
	
	public static GraphStyle getDefaultGraphStyle(final int funcIndex, final boolean asLine, DrawingSupplier drawingSupplier) {
		if(asLine) {
			lastShape = funcIndex % 5 == 0? drawingSupplier.getNextShape() : lastShape;
			return new GraphStyle(lastShape, new BasicStroke(2.0f), getDefaultColor(funcIndex));
		} else {
			return new GraphStyle(drawingSupplier.getNextShape(), null, getDefaultColor(funcIndex));
		}
	}
	
	public static Color getColor(final double red, final double green, final double blue) {
		return new java.awt.Color((int)(red*255), (int)(green*255), (int)(blue*255));
	}
	
	public static Color getColor(final int red, final int green, final int blue) {
		return new java.awt.Color(red, green, blue);
	}

	// Super methods:
	@Override
	public EasyPlot2D setTitle(String title) {
		this.title = title;
		super.setTitle(title);
		return this;
	}

	@Override
	public EasyPlot2D setXAxisLabel(String xAxisLabel) {
		xAxisName = xAxisLabel;
		super.setXAxisLabel(xAxisLabel);
		return this;
	}

	@Override
	public EasyPlot2D setYAxisLabel(String yAxisLabel) {
		yAxisName = yAxisLabel;
		super.setYAxisLabel(yAxisLabel);
		return this;
	}

	@Override
	public EasyPlot2D setXRange(double xmin, double xmax) {
		xMin = xmin;
		xMax = xmax;
		super.setXRange(xmin, xmax);
		return this;
	}

	@Override
	public EasyPlot2D setYRange(double ymin, double ymax) {
		yMin = ymin;
		yMax = ymax;
		super.setYRange(ymin, ymax);
		return this;
	}

	@Override
	public EasyPlot2D setIsLegendVisible(Boolean isLegendVisible) {
		legend = isLegendVisible;
		super.setIsLegendVisible(isLegendVisible);
		return this;
	}

	private static class Plotable2DCloneable implements Plotable2D {

		final String name;
		final List<Point2D> series;
		final GraphStyle style;
		final NumberAxis domainAxis;
		final NumberAxis rangeAxis;

		final Kind kind;
		public static enum Kind {
			FUNCTION,
			POINTS,
			THIS,
			UNKNOWN
		}

		private Plotable2DCloneable(String name, List<Point2D> series, GraphStyle style, NumberAxis domainAxis, NumberAxis rangeAxis, Kind kind) {
			this.name = name;
			this.series = series;
			this.style = style;
			this.domainAxis = domainAxis;
			this.rangeAxis = rangeAxis;
			this.kind = kind;
		}

		public Plotable2DCloneable(String name, List<Point2D> series, GraphStyle style, NumberAxis domainAxis, NumberAxis rangeAxis) {
			this(name, series, style, domainAxis, rangeAxis, Kind.THIS);
		}

		public Plotable2DCloneable(String name, List<Point2D> series, GraphStyle style) {
			this(name, series, style, null, null);
		}
		Plotable2DCloneable(Plotable2D plotable) {
			this.name = plotable.getName();
			this.series = plotable.getSeries();
			this.style = plotable.getStyle();
			this.domainAxis = plotable.getDomainAxis();
			this.rangeAxis = plotable.getRangeAxis();
			this.kind = getKind(plotable);
		}

		public static Plotable2DCloneable of(Plotable2D plotable) {
			return new Plotable2DCloneable(plotable);
		}

		public static Plotable2DCloneable of(String name, List<Point2D> series, GraphStyle style, NumberAxis domainAxis, NumberAxis rangeAxis, Plotable2DCloneable.Kind kind) {
			return new Plotable2DCloneable(name, series, style, domainAxis, rangeAxis, kind);
		}

		public static List<Plotable2D> listOf(final List<Plotable2D> plotables, DrawingSupplier leDrawingSupplier) {
			Function<Integer, Plotable2D> mapper = (i) -> {
				if(plotables.get(i).getStyle() == null) {
					return of(plotables.get(i).getName(), plotables.get(i).getSeries(), getDefaultGraphStyle(i, true, leDrawingSupplier), plotables.get(i).getDomainAxis(), plotables.get(i).getRangeAxis(), getKind(plotables.get(i)));
				}
				else
					return of(plotables.get(i));
			};
			return IntStream.range(0, plotables.size()).mapToObj(mapper::apply).collect(Collectors.toList());
		}

		public Plotable2DCloneable clone() {
			return new Plotable2DCloneable(this);
		}

		public Plotable2DCloneable getCloneWithModifiedName(String newName) {
			return new Plotable2DCloneable(newName, getSeries(), getStyle(), getDomainAxis(), getRangeAxis(), getKind());
		}

		public Plotable2DCloneable getCloneWithModifiedSeries(List<Point2D> newSeries) {
			return new Plotable2DCloneable(getName(), newSeries, getStyle(), getDomainAxis(), getRangeAxis(), getKind());
		}

		public Plotable2DCloneable getCloneWithModifiedStyle(GraphStyle newStyle) {
			return new Plotable2DCloneable(getName(), getSeries(), newStyle, getDomainAxis(), getRangeAxis(), getKind());
		}

		public Plotable2DCloneable getCloneWithModifiedDomainAxis(NumberAxis newDomainAxis) {
			return new Plotable2DCloneable(getName(), getSeries(), getStyle(), newDomainAxis, getRangeAxis(), getKind());
		}

		public Plotable2DCloneable getCloneWithModifiedRangeAxis(NumberAxis newRangeAxis) {
			return new Plotable2DCloneable(getName(), getSeries(), getStyle(), getDomainAxis(), newRangeAxis, getKind());
		}

		@Override
		public List<Point2D> getSeries() {
			return series;
		}

		@Override
		public GraphStyle getStyle() {
			return style;
		}

		@Override
		public NumberAxis getDomainAxis() {
			return domainAxis;
		}

		@Override
		public NumberAxis getRangeAxis() {
			return rangeAxis;
		}

		@Override
		public String getName() {
			return name;
		}

		public Kind getKind() {
			return kind;
		}

		private static Kind getKind(Plotable2D origin) {
			if(origin instanceof PlotableFunction2D)
				return Kind.FUNCTION;
			else if(origin instanceof PlotablePoints2D)
				return Kind.POINTS;
			else if(origin instanceof Plotable2DCloneable pl)
				return pl.kind;
			else
				return Kind.UNKNOWN;

		}

	}
}
