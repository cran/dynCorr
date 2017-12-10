# dynCorr 1.1.0
* Expanded functionality to incorporate dynamic correlation estimation for datasets with non-integer valued indep variable. 
* A column is added to return value of bootstrapCI to give direct visual comparison of cor estimation and confidence interval estimation. 
* Table data.summary is added to return values of both bootstrapCI and dynamicCorrelation to summarize heterogeneity of dataset. Three columns included in this table are: sample.size (number of observations used in the calculation), min.max.time (max common observation time), max.max.time (max individual observation time in the sample).
* Three new variables are introduced to both dynamicCorrelation and bootstrapCI: points.length, points.by, and min.obs. See help page for more detail.


# dynCorr 1.0.0

* Added a `NEWS.md` file to track changes to the package.
* Reformated code, and removed unnecessary code (e.g., usage of attach and detach), 
    to adhere to modern practice/standards.
* Fixed a bug in dynamicCorrelation (was applicable only to input with byOrder specified, lag < 0, and  by.deriv.only = FALSE).
* Fixed a typo regarding parentheses usage in bootstrapCI.



