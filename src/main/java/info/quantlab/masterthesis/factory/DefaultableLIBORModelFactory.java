package info.quantlab.masterthesis.factory;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.function.BiFunction;

import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceModel;
import info.quantlab.masterthesis.defaultablecovariancemodels.DefaultableLIBORCovarianceWithGuaranteedPositiveSpread;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORFromSpreadDynamic;
import info.quantlab.masterthesis.defaultablelibormodels.DefaultableLIBORMarketModel;
import info.quantlab.masterthesis.multilibormodels.MultiLIBORVectorModel;
import info.quantlab.masterthesis.process.MilsteinSchemeAnalyticDerivative;
import info.quantlab.masterthesis.process.MilsteinSchemeFiniteDifference;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.CalibrationProduct;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.BlendedLocalVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.DisplacedLocalVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.model.ProcessModel;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.opencl.montecarlo.RandomVariableOpenCLFactory;
import net.finmath.randomnumbers.MersenneTwister;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import net.finmath.util.TriFunction;

/**
 * This class is used for a flexible construction of a DefaultableLIBORMarketModel. All model specific values (i.e. fields) have default
 * values (except for values unused in the default state) and can be modified through the {@link #setProperties} method.
 * 
 * @author Markus Parhofer
 * @version 0.9
 */
public class DefaultableLIBORModelFactory {

	/**
	 * Model implementation used for the non defaultable covariance structure.
	 * <br>
	 * Values:  <br>
	 * <b>NORMAL</b>, <br>
	 * <b>DISPLACED</b>, <br>
	 * <b>BLENDED</b>
	 * @author Markus Parhofer
	 * @version 1.0
	 */
	public enum CovarianceModel {
		NORMAL,
		DISPLACED,
		BLENDED
	}
	
	/**
	 * Scheme Implementation used for the MonteCarloProcess.
	 * <br>
	 * Values:  <br>
	 * <b>EULER</b>, <br>
	 * <b>EULER_FUNCTIONAL</b>, <br>
	 * <b>MILSTEIN_ANALYTIC</b>, <br>
	 * <b>MILSTEIN_FDFORWARD</b>, <br>
	 * <b>MILSTEIN_FDBACKWARD</b>, <br>
	 * <b>MILSTEIN_FDCENTRAL</b>
	 * @author Markus Parhofer
	 * @version 1.0
	 */
	public enum Scheme {
		EULER,
		EULER_FUNCTIONAL,
		MILSTEIN_ANALYTIC,
		MILSTEIN_FDFORWARD,
		MILSTEIN_FDBACKWARD,
		MILSTEIN_FDCENTRAL,
	}
	
	
	// General Parameters:
	private double[] fixingTimes = { 0.5, 1.0, 2.0, 4.0, 8.0, 25.0 };
	private double liborPeriodLength = 2.0;
	private int numberOfLiborPeriods = 5;
	private String stateSpace = "NORMAL";
	private String measure = "SPOT";

	
	// Non Defaultable Model Parameters:
	private CovarianceModel covarianceModel = CovarianceModel.NORMAL;
	private int numberOfFactors = 5; // No Factor reduction
	private double[] volatilityParams = { 0.1, 0.0, 0.25, 0.1 };
	private double correlationDecayParam = 0.2;
	private double displacement = 0.0;
	private double[] initialRatesNonDefaultable = { 0.035, 0.043, 0.05, 0.041, 0.035, 0.02 };
	
	
	// Defaultable Model Parameters:
	private int numberOfExtraFactors = 2;
	private double[] initialRatesDefaultable = { 0.04, 0.048, 0.055, 0.046, 0.04, 0.025 };
	private String simulationModel;
	private String stateSpaceOfSpread;
	private int freeParamsSeed = 111;
	private double freeParamsRange = 0.5;
	private BiFunction<Integer, Double, double[][]> freeParamsGenerator = (seed, range) -> {
		return generateFreeParamsDefault(seed, range, numberOfLiborPeriods, numberOfExtraFactors);
	};
	
	
	// Simulation parameters:
	private double simulationTimeDelta = 0.001;
	private int numberOfPaths = 1000;
	private int brownianMotionSeed = 7824;
	private Scheme numericalScheme = Scheme.EULER;
	private double finiteDifferenceDelta = 1E-6;
	private TriFunction<Integer, Double, RandomVariable[], RandomVariable[]> analyticDifferentialFactorLoadings = null;
	
	
	public DefaultableLIBORModelFactory() {}

	public Map<String, Object> getProperties() {
		Map<String, Object> properties = new HashMap<>();
		
		
		properties.put("fixingTimes", fixingTimes);
		properties.put("liborPeriodLength", liborPeriodLength);
		properties.put("numberOfLiborPeriods", numberOfLiborPeriods);
		properties.put("stateSpace", stateSpace);
		properties.put("measure", measure);
		
		
		properties.put("covarianceModel", covarianceModel);
		properties.put("numberOfFactors", numberOfFactors);
		properties.put("volatilityParams", volatilityParams);
		properties.put("correlationDecayParam", correlationDecayParam);
		properties.put("displacement", displacement);
		properties.put("initialRatesNonDefaultable", initialRatesNonDefaultable);
		
		
		properties.put("numberOfExtraFactors", numberOfExtraFactors);
		properties.put("initialRatesDefaultable", initialRatesDefaultable);
		properties.put("simulationModel", simulationModel);
		properties.put("stateSpaceOfSpread", stateSpaceOfSpread);
		properties.put("freeParamsSeed", freeParamsSeed);
		properties.put("freeParamsRange", freeParamsRange);
		properties.put("freeParamsGenerator", freeParamsGenerator);
		
		
		properties.put("simulationTimeDelta", simulationTimeDelta);
		properties.put("numberOfPaths", numberOfPaths);
		properties.put("brownianMotionSeed", brownianMotionSeed);
		properties.put("numericalScheme", numericalScheme);
		properties.put("finiteDifferenceDelta", finiteDifferenceDelta);
		properties.put("analyticDifferentialFactorLoadings", analyticDifferentialFactorLoadings);
		
		
		return properties;
	}
	
	/**
	 * Valid Keys for the properties are:
	 * <ul>
	 * <li><b>fixingTimes</b> (<code>double[]</code>): Fixing times associated with the initial forward rates.</li>
	 * <li><b>liborPeriodLength</b> (<code>double</code>): length of the LIBOR periods.</li>
	 * <li><b>numberOfLiborPeriods</b> (<code>int</code>): number of LIBOR periods.</li>
	 * <li><b>stateSpace</b> (<code>String</code>): the statespace of the defaultable and non defaultable model.</li>
	 * <li><b>measure</b> (<code>String</code>): the measure of the defaultable and non defaultable model.</li>
	 * 
	 * <li><b>covarianceModel</b> (<code>String</code> or {@link DefaultableLIBORModelFactory.CovarianceModel}): covariance model implementation used for the non defaultable model.</li>
	 * <li><b>numberOfFactors</b> (<code>int</code>): number of factors for the non defaultable model.</li>
	 * <li><b>volatilityParams</b> (<code>double[]</code>): parameters for the volatility structure of the non defaultable model.</li>
	 * <li><b>correlationDecayParam</b> (<code>double</code>): correlation decay parameter for the non defaultable correlation model.</li>
	 * <li><b>displacement</b> (<code>double</code>): displacement parameter for the displaced or blended non defaultable covariance model.</li>
	 * <li><b>initialRatesNonDefaultable</b> (<code>double[]</code>): initial forward rates of the non defaultable model associated with fixingTimes.</li>
	 * 
	 * <li><b>numberOfExtraFactors</b> (<code>int</code>): number of additional factors for the covariance structure of the defaultable model.</li>
	 * <li><b>initialRatesDefaultable</b> (<code>double[]</code>): initial forward rates of the defaultable model associated with fixingTimes.</li>
	 * <li><b>simulationModel</b> (<code>String</code>): name of the value to model by the MonteCarloProcess to get the defaultable LIBOR rates (i.e. SPREADS or LIBORS).</li>
	 * <li><b>stateSpaceOfSpread</b> (<code>String</code>): the statespace of the spreads.</li>
	 * <li><b>freeParamsSeed</b> (<code>int</code>): the seed for the free parameter matrix generator.</li>
	 * <li><b>freeParamsRange</b> (<code>double</code>): the range for the free parameter matrix (all free paramters will be in [-range, range]).</li>
	 * <li><b>freeParamsGenerator</b> (<code>BiFunction{@literal <}Integer, Double, double[][]{@literal >}</code>): the function used to generate the free parameter matrix. Input values are (freeParamsSeed, freeParamsRange).</li>
	 * 
	 * <li><b>simulationTimeDelta</b> (<code>double</code>): the discretization time step that is used for simulation and the construction of the covariance structures.</li>
	 * <li><b>numberOfPaths</b> (<code>int</code>): number of paths used for for the simulation.</li>
	 * <li><b>brownianMotionSeed</b> (<code>int</code>): seed used for the construction of the BrownianMotionFromMersenneRandomNumbers.</li>
	 * <li><b>numericalScheme</b> (<code>String</code> or {@link DefaultableLIBORModelFactory.Scheme}): scheme used to simulate the model.</li>
	 * <li><b>finiteDifferenceDelta</b> (<code>double</code>): the delta used to estimate the differential of the factor loadings in a finite difference Milstein scheme.</li>
	 * <li><b>analyticDifferentialFactorLoadings</b> (<code>TriFunction{@literal <}Integer, Double, RandomVariable[], RandomVariable[]{@literal >}</code>): function used for the differential of the factor loadings. Input values are (componentIndex, time, LIBORrealizations).</li>
	 * </ul>
	 * @apiNote All properties are unchecked. If two properties are incompatible this method will not throw an exception.
	 * @param properties Map with specified properties. see above for options
	 */
	@SuppressWarnings("unchecked")
	public void setProperties(Map<String, Object> properties) {
		if(properties != null) {
			Set<String> keys = properties.keySet();
			for(String key: keys) {
				switch(key.toLowerCase()) {
				case "fixingtimes":
					fixingTimes = (double[]) properties.get(key);
					break;
				case "liborperiodlength":
					liborPeriodLength = (double) properties.get(key);
					break;
				case "numberofliborperiods":
					numberOfLiborPeriods = (int) properties.get(key);
					break;
				case "statespace":
					stateSpace = (String) properties.get(key);
					break;
				case "measure":
					measure = (String) properties.get(key);
					break;
					
				case "covariancemodel":
					if(properties.get(key) instanceof CovarianceModel covModel) {
						covarianceModel = covModel;
					}
					else if(properties.get(key) instanceof String stringCovModel) {
						covarianceModel = CovarianceModel.valueOf(stringCovModel);
					}
					break;
				case "numberoffactors":
					numberOfFactors = (int) properties.get(key);
					break;
				case "volatilityparams":
					volatilityParams = (double[]) properties.get(key);
					break;
				case "correlationdecayparam":
					correlationDecayParam = (double) properties.get(key);
					break;
				case "displacement":
					displacement = (double) properties.get(key);
					break;
				case "initialratesnondefaultable":
					initialRatesNonDefaultable = (double[]) properties.get(key);
					break;
					
				case "numberofextrafactors":
					numberOfExtraFactors = (int) properties.get(key);
					break;
				case "initialratesdefaultable":
					initialRatesDefaultable = (double[]) properties.get(key);
					break;
				case "simulationmodel":
					simulationModel = (String) properties.get(key);
					break;
				case "statespaceofspread":
					stateSpaceOfSpread = (String) properties.get(key);
					break;
				case "freeparamsseed":
					freeParamsSeed = (int) properties.get(key);
					break;
				case "freeparamsrange":
					freeParamsRange = (double) properties.get(key);
					break;
				case "freeparamsgenerator":
					freeParamsGenerator = (BiFunction<Integer, Double, double[][]>) properties.get(key);
					break;
					
				case "simulationtimedelta":
					simulationTimeDelta = (double) properties.get(key);
					break;
				case "numberofpaths":
					numberOfPaths = (int) properties.get(key);
					break;
				case "brownianmotionseed":
					brownianMotionSeed = (int) properties.get(key);
					break;
				case "numericalscheme":
					if(properties.get(key) instanceof Scheme scheme) {
						numericalScheme = scheme;
					}
					else if(properties.get(key) instanceof String stringScheme) {
						numericalScheme = Scheme.valueOf(stringScheme);
					}
					break;
				case "finitedifferencedelta":
					finiteDifferenceDelta = (double) properties.get(key);
					break;
				case "analyticdifferentialfactorloadings":
					analyticDifferentialFactorLoadings = (TriFunction<Integer, Double, RandomVariable[], RandomVariable[]>) properties.get(key);
					break;
				}
			}
		}

	}
	
	
	public LIBORMarketModel createBaseModel() throws CalculationException {
		// Set LIBOR times
		TimeDiscretization liborPeriods = new TimeDiscretizationFromArray(0.0, numberOfLiborPeriods, liborPeriodLength); // Fixing time
		
		// Set initial forward curves
		final ForwardCurve nonDefaultableForwards = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"nonDefaultableForwardCurve",
				fixingTimes,
				initialRatesNonDefaultable,
				liborPeriodLength);
		
		/*
		 * Create a simulation time discretization. We save space by just simulationg to the last fixing time.
		 */
		final int numberOfTimes = (int)(liborPeriods.getTime(liborPeriods.getNumberOfTimeSteps() - 1)/simulationTimeDelta);
		final TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(0.0, numberOfTimes, simulationTimeDelta);

		/*
		 * Create volatility structure v[i][j] = sigma_j(t_i) and correlation model rho_{i,j} = exp(-a * abs(T_i-T_j))
		 */
		
		final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(timeDiscretization, liborPeriods, volatilityParams[0], volatilityParams[1], volatilityParams[2], volatilityParams[3], true);
		final LIBORCorrelationModelExponentialDecay correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretization, liborPeriods, numberOfFactors,	correlationDecayParam, true);
		
		AbstractLIBORCovarianceModelParametric baseCovarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization, liborPeriods, volatilityModel, correlationModel);
		switch(covarianceModel) {
		case NORMAL:
			break;
		case BLENDED:
			baseCovarianceModel = new BlendedLocalVolatilityModel(baseCovarianceModel, null, displacement, true);
			break;
		case DISPLACED:
			baseCovarianceModel = new DisplacedLocalVolatilityModel(baseCovarianceModel, displacement, true);
			break;
		}
		
		final Map<String, String> properties = new HashMap<>();
		properties.put("measure", measure);
		properties.put("stateSpace", stateSpace);

		// Empty array of calibration items - hence, model will use given covariance
		final CalibrationProduct[] calibrationItems = new CalibrationProduct[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */
		return LIBORMarketModelFromCovarianceModel.of(
				liborPeriods, 		/* LIBORPeriodDiscretization */
				null,			 	/* analyticModel */
				nonDefaultableForwards, /* ForwardCurve */
				new DiscountCurveFromForwardCurve(nonDefaultableForwards),  /* DiscountCurve */
				new RandomVariableFromArrayFactory(true),  /* RV Factory */
				baseCovarianceModel,  /* covarianceModel */
				calibrationItems,  /* calibrationItems */
				properties);
	}

	public DefaultableLIBORMarketModel createDefaultableModel(LIBORMarketModel baseModel) {
		final ForwardCurve defaultableForwards = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"defaultableForwardCurve",
				fixingTimes,
				initialRatesDefaultable,
				liborPeriodLength);
		
		final Map<String, String> properties = new HashMap<>();
		properties.put("measure", measure);
		properties.put("stateSpace", stateSpace);
		properties.put("simulationModel", simulationModel);
		properties.put("stateSpaceOfSpread", stateSpaceOfSpread);
		
		final DefaultableLIBORCovarianceModel defaultableCovariance = new DefaultableLIBORCovarianceWithGuaranteedPositiveSpread(baseModel.getCovarianceModel(), freeParamsGenerator.apply(freeParamsSeed, freeParamsRange));
		//final DefaultableLIBORCovarianceModel defaultableCovariance = new DefaultableLIBORCovarianceWithInitialUndefaultableCovariance(baseCovarianceModel, defaultableForwards, nonDefaultableForwards, numberOfExtraFactors, 0.5);
				
		return new DefaultableLIBORFromSpreadDynamic(baseModel, defaultableCovariance, defaultableForwards, properties);
	}
	
	
	public DefaultableLIBORMarketModel createDefaultableModel() throws CalculationException {
		return createDefaultableModel(createBaseModel());
	}
	
	
	public MultiLIBORVectorModel createMultiModel(LIBORMarketModel baseModel, Map<String, Object>[] defaultableProperties) {
		DefaultableLIBORMarketModel[] defModels = new DefaultableLIBORMarketModel[defaultableProperties.length];
		for(int modelIndex = 0; modelIndex < defaultableProperties.length; modelIndex++) {
			 DefaultableLIBORModelFactory defModelFactory = clone();
			 defModelFactory.setProperties(defaultableProperties[modelIndex]);
			 defModels[modelIndex] = defModelFactory.createDefaultableModel(baseModel);
		}
		return new MultiLIBORVectorModel(defModels, baseModel);
	}
	
	
	public MonteCarloProcess createNumericalScheme(ProcessModel model) throws NullPointerException {
		TimeDiscretization timeDiscretization;
		{
			var objectModel = (Object)model;
			if(objectModel instanceof LIBORMarketModel liborModel) {
				timeDiscretization = liborModel.getCovarianceModel().getTimeDiscretization();
			}
			else {
				final int numberOfTimes = (int)((double)(numberOfLiborPeriods - 1) * liborPeriodLength/simulationTimeDelta);
				timeDiscretization = new TimeDiscretizationFromArray(0.0, numberOfTimes, simulationTimeDelta);
			}
		}
		
		final BrownianMotion brownianMotion = new BrownianMotionFromMersenneRandomNumbers(timeDiscretization, model.getNumberOfFactors(), numberOfPaths, brownianMotionSeed);

        return switch (numericalScheme) {
            case EULER ->
                    new EulerSchemeFromProcessModel(model, brownianMotion, EulerSchemeFromProcessModel.Scheme.EULER);
            case EULER_FUNCTIONAL ->
                    new EulerSchemeFromProcessModel(model, brownianMotion, EulerSchemeFromProcessModel.Scheme.EULER_FUNCTIONAL);
            case MILSTEIN_ANALYTIC -> {
                if (analyticDifferentialFactorLoadings == null)
                    throw new NullPointerException("For the Milstein Scheme with analytic differential the analytic differential of the Factor Loadings must be set!");
                yield new MilsteinSchemeAnalyticDerivative(model, brownianMotion, analyticDifferentialFactorLoadings);
            }
            case MILSTEIN_FDBACKWARD -> {
                MilsteinSchemeFiniteDifference.FiniteDifferenceMethod methodBack = MilsteinSchemeFiniteDifference.FiniteDifferenceMethod.BACKWARD;
                yield new MilsteinSchemeFiniteDifference(model, brownianMotion, methodBack, finiteDifferenceDelta);
            }
            case MILSTEIN_FDFORWARD -> {
                MilsteinSchemeFiniteDifference.FiniteDifferenceMethod methodFor = MilsteinSchemeFiniteDifference.FiniteDifferenceMethod.FORWARD;
                yield new MilsteinSchemeFiniteDifference(model, brownianMotion, methodFor, finiteDifferenceDelta);
            }
            case MILSTEIN_FDCENTRAL -> {
                MilsteinSchemeFiniteDifference.FiniteDifferenceMethod methodCent = MilsteinSchemeFiniteDifference.FiniteDifferenceMethod.CENTRAL;
                yield new MilsteinSchemeFiniteDifference(model, brownianMotion, methodCent, finiteDifferenceDelta);
            }
        };
	}
	
	public MonteCarloProcess createNumericalScheme() throws NullPointerException, CalculationException {
		return createNumericalScheme(createDefaultableModel());
	}
	
	public DefaultableLIBORModelFactory clone() {
		DefaultableLIBORModelFactory newFactory = new DefaultableLIBORModelFactory();
		newFactory.setProperties(getDefaultProperties());
		return newFactory;
	}
	
	public static Map<String, Object> getDefaultProperties() {
		DefaultableLIBORModelFactory newFactory = new DefaultableLIBORModelFactory();
		return newFactory.getProperties();
	}
	
	public static double[][] generateFreeParamsDefault(final int seed, final double range, final int liborPeriods, final int extraFactors) {
		double[][] matrix = new double[liborPeriods][extraFactors];
		MersenneTwister randomGenerator = new MersenneTwister(seed);
		for (int row = 0; row < matrix.length; row++) {
			for (int col = 0; col < matrix[row].length; col++) {
				matrix[row][col] = randomGenerator.nextDoubleFast() * range * 2.0d - range;
			}
		}
		return matrix;
	}
}
