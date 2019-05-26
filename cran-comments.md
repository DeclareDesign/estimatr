## Test environments
* local OS X install
* travis-ci ubuntu (3.4, 3.5, 3.6)
* travis-ci OS X (3.4, 3.5)
* Appveyor windows (3.4, 3.5)

Passes with only notes about package size.

However, the main purpose of this update is to respond to CRAN errors due to updates to clubSandwich. These have been resolved.

Unfortunately, other outstanding errors with some memory issues on the CRAN solaris boxes have not been resolved due to an inability to test package updates in a stable Solaris environment (r-hub's facilities are not working). We are still working on this issue.
