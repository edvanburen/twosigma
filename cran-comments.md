**COMMENTS TO CRAN**

This is a resubmission of the twosigma package, version 1.0.2. Based on the feedback I received from CRAN, I made the following changes to my submission:

1. Shortened the length of the title of the package

2. Added two references relevant to the package inside the description field of DESCRIPTION

3. Added examples inside of all eight exported functions

4. Added tests for all eight exported functions in "tests/testthat/test-twosigma.R"

5. Removed the case in which a user's global options were changed from the function twosigmag.R
(commented line 99), and removed the function twosigmag_ttest.R entirely from the package (also had previously modified a user's global options)
