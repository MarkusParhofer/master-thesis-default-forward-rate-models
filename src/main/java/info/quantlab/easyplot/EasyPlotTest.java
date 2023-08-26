package info.quantlab.easyplot;

import java.awt.BasicStroke;
import java.io.IOException;
import java.util.function.DoubleUnaryOperator;

import net.finmath.plots.GraphStyle;
import net.finmath.plots.Named;

public class EasyPlotTest {

	public static void main(String[] args) {
		
		Named<DoubleUnaryOperator> funcs0 = new Named<DoubleUnaryOperator>("Sinus", new DoubleUnaryOperator() {
			@Override
			public double applyAsDouble(double operand) {
				return Math.sin(operand);
			}
		});
		
		Named<DoubleUnaryOperator> funcs1 = new Named<DoubleUnaryOperator>("Cosinus", new DoubleUnaryOperator() {
			@Override
			public double applyAsDouble(double operand) {
				return Math.cos(operand);
			}
		});
		
		
		Named<DoubleUnaryOperator> funcs2 = new Named<DoubleUnaryOperator>("Absolut Wert", new DoubleUnaryOperator() {
			@Override
			public double applyAsDouble(double operand) {
				return Math.abs(operand);
			}
		});
		
		EasyPlot2D myPlot = new EasyPlot2D(0.0, 5.0, 300, funcs0, new GraphStyle(null, new BasicStroke(2.0f), EasyPlot2D.getColor(0.5,  0.5, 1.0)));
		
		myPlot.addPlot(0.0, 4.0, 300, funcs1);
		myPlot.show();
		
		try {
			System.in.read();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		myPlot.addPlot(0.0, 5.0, 300, funcs2);
		myPlot.changePlotColor(1, EasyPlot2D.getColor(1.0, 0.0, 0.0));
		myPlot.setIsLegendVisible(true);
	}

}
