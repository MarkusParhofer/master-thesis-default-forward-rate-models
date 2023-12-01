
# master-thesis-default-forward-rate-models

This is an OOP implementation for the construction of defaultable forward rate models and valuation examples through these.
The models follow the theory of the Paper [Defaultable Discrete Forward Rate Model with Covariance Structure Guaranteeing Positive Credit Spreads](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3667878) authored by Professor Dr. Christian Fries.
The implementation in java and embedding into the Finmath Library of Prof. Fries is part of the Master's thesis "Defaultable forward rate models for the valuation of loans with behavioral aspects" authored by Markus Parhofer which will soon be available here.


**Author:** Markus Parhofer

**Supervision:** Prof. Dr. Christian Fries, Prof. Dr. Andrea Mazzon

**Version: 0.0.1 (currently in development)**

## Important Note
The models are ***currently in Development*** and not ready for any kind of real world application. Use at your own risk.
Generally no warranty for the correctness of the implementation is given. For further details read the [license](./LICENSE.txt) and [notice](./NOTICE.txt).


## Requirements
- Java SE17
- Finmath Library
- Finmath Plot Library

## Structure

The classes are embedded into the framework of the Finmath library of Prof. Dr. Christian Fries. 
### Packages:
- **[info.quantlab.masterthesis.defaultablelibormodels](./src/main/java/info/quantlab/masterthesis/defaultablelibormodels):** This is the heart of the repository. Here we have the classes that implement the defaultable LIBOR Market models. They can be used as a plug in for the <code>EulerSchemeFromProcessModel</code> class of the finmath lib. I recommend using <code>DefaultableLIBORFromSpreadDynamic</code>.
- **[info.quantlab.masterthesis.defaultablecovariancemodels](./src/main/java/info/quantlab/masterthesis/defaultablecovariancemodels):** This package holds the covariance models for plugging into the <code>DefaultableLiborMarketModel</code>. I recommend using <code>DefaultableLIBORCovarianceWithGuaranteedPositiveSpread</code>. Other models have not been tested or are not working correctly.
- **[info.quantlab.masterthesis.multilibormodels](./src/main/java/info/quantlab/masterthesis/multilibormodels):** Holds the class for the usage with multiple defaultable libor models.
- **[info.quantlab.masterthesis.process](./src/main/java/info/quantlab/masterthesis/process):** Implements a Milstein-Scheme for the simulation of solutions to stochastic differential equations given as <code>ProcessModels</code>.
- **[info.quantlab.masterthesis.products](./src/main/java/info/quantlab/masterthesis/products):** Implements some defaultable products that are designed for the usage with the <code>DefaultableLiborMarketModel</code> and the <code>MultiLiborVectorModel</code>
- **[info.quantlab.masterthesis.functional](./src/main/java/info/quantlab/masterthesis/functional):** Classes with static functions for flexible use of <code>BrownianMotion</code> and <code>MonteCarloProcess</code> as well as <code>RandomVariable</code>.
- **[info.quantlab.easyplot](./src/main/java/info/quantlab/easyplot):** Has a class to plot various content with easy access to the plots styles.
- **[info.quantlab.debug](src/main/java/info/quantlab/debug):** Has tools and classes that are good for debugging.
- **[info.quantlab.masterthesis.legacy](./src/main/java/info/quantlab/masterthesis/legacy):** In here are old models that proved not to be very usable, but might be helpful for later restructuring.

### Tests:
The only meaningful test file at the moment is [ModelFromSpreadTest](./src/test/java/info/quantlab/masterthesis/defaultablemodels/testing/ModelFromSpreadTest.java).
#### Setting Model Parameters:
- **Automatic Config**: To compare the performance on different configurations one can specify some main parameters through the constructor values.
- **Manual Config setting**: To have a more advanced control over config specifications one can change the values in the map that is constructed in the method <code>getValuationModel(...)</code>.

#### Plots:
The test automatically constructs plots.
- **Saving the plots**: To save the plots one can use the boolean <code>savePlots</code> (default: <code>true</code>). This will generate a folder in [Graphs](./Graphs) with the time and date stamp and the title of the run configurations.
- **Showing the plots**: The plots will only be shown if <code>savePlots</code> is set to <code>false</code>. Only the plot with the averages will be shown. Note that in junit tests the plots automatically close once the tests are completed.