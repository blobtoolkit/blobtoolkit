# Template Python3 project repository configured for Continuous Integration

```
git clone --bare https://github.com/rjchallis/test.git
```

```
cd test.git
git push --mirror https://github.com/exampleuser/new-repository.git
```

```
cd ..
rm -rf test.git
```

[![Build Status](https://travis-ci.org/rjchallis/test.svg?branch=master)](https://travis-ci.org/rjchallis/test)
[![Coverage Status](https://coveralls.io/repos/github/rjchallis/test/badge.svg?branch=master)](https://coveralls.io/github/rjchallis/test?branch=master)
