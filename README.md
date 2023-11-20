# master-thesis-default-forward-rate-models
Master's thesis to implement defaultable forward rate models for the valuation of loans with behavioral aspects.

**Author:** Markus Parhofer

**Supervision:** Prof. Dr. Christian Fries, Prof. Dr. Andrea Mazzon

## Requirements
- Java SE17
- Finmath Library
- Finmath Plot Library

## Structure

The classes are embedded into the framework of the Finmath library of Prof. Dr. Christian Fries. 
Packages:
- **info.quantlab.masterthesis.defaultablelibormodels:** This is the heart of the repository. Here we have the classes that implement the defaultable LIBOR Market models. They can be used as a plug in for the EulerSchemeFromProcessModel class of the finmath lib. I recommend using "DefaultableLIBORFromSpreadDynamic".
- **info.quantlab.masterthesis.defaultablecovariancemodels:** This package holds the covariance models for plugging into the DefaultableLiborMarketModel. I recommend using "DefaultableLIBORCovarianceWithGuaranteedPositiveSpread". Other models have not been tested or are not working correctly.
- **info.quantlab.masterthesis.products:** 
- **info.quantlab.masterthesis.functional:** Class with static functions for flexible use of BrownianMotions and MonteCarloProcesses.
- **info.quantlab.easyplot:** Has Class to plot various content with easy access to the plots styles.
- **info.quantlab.debug:** Has tools and classes that are good for debugging.
- **info.quantlab.masterthesis.legacy:** In here are old models that proofed to not be very useable, but might be helpful for later restructuring.
