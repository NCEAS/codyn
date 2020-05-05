# Contributing to codyn

:tada: First off, thanks for contributing!:tada:

- [Types of contributions](#types-of-contributions)
- [Pull Requests](#pull-requests)
- [Development Workflow](#development-workflow)
- [Release process](#release-process)
- [Testing](#testing)
- [Code style](#code-style)
- [Contributor license agreement](#contributor-license-agreement)

## Types of contributions

We welcome all types of contributions, including bug fixes, feature enhancements,
bug reports, documentation, graphics, and many others.  You might consider contributing by:

- Report a bug or request a new feature in our [issue tracker](https://github.com/DataONEorg/scythe/issues)
- Fix a bug and contribute the code with a Pull Request
- Write or edit some documentation
- Develop a screencast tutorial
- Answer questions on our mailing list or Slack team
- ...

codyn is an open source project, and we welcome full
participation in the project.  Contributions are reviewed and suggestions are
made to increase the value of codyn to the community.  We strive to
incorporate code, documentation, and other useful contributions quickly and
efficiently while maintaining a high-quality product.

## Pull Requests
We use the pull-request model for contributions. See [GitHub's help on pull-requests](https://help.github.com/articles/about-pull-requests/).

In short:

- add an [issue](https://github.com/DataONEorg/scythe/issues) describing your planned changes, or add a comment to an existing issue;
- on GitHub, fork the [codyn repository](https://github.com/NCEAS/codyn)
- on your computer, clone your forked copy of the repository
- checkout the `develop` branch and commit your changes
- push your branch to your forked repository, and submit a pull-request
- go through the review process, making changes until your pull-request can be merged
- describe your changes in the issue, and close your issue.

## Development Workflow

Development is managed through the git repository at https://github.com/NCEAS/codyn.  The repository is organized into several branches, each with a specific purpose. In GitHub, users will normally see the `master` branch, which therefore should reflect the current stable release of the package (rather than confusing people with in-progress proposed changes that are not yet released). Consequently, we will use a [GitFlow](https://nvie.com/posts/a-successful-git-branching-model/)-inspired release model in which the `master` branch always reflects the current stable release, a `develop` branch is used for merging finished proposals being prepared for the next release, and `feature` branches are used for creating changes to implement specific enhancements and bugfixes.  For minor changes that don't require review, such as spelling corrections, grammatical rewording, etc., maintainers can commit changes directly to the `develop` branch, and other contributors can do a pull request directly against the `develop` branch.  The use of feature branches is really focused on managing proposals that need discussion, review, and a decision by the team. Maintainers will make judgement calls on whether an feature branch is needed, and might convert contributed pull requests to a feature branch if they determine that one  is needed.

**master**. The `master` branch represents the stable branch that contains the current release. Anybody installing the package will, therefore, by default get the stable release.  

**develop**. Development takes place on the `develop` branch, which represents the 
changes being planned for the next release. This is where the next release is fully tested,
and all feature changes should be merged onto `develop` and tested.


When a set of features are mature and tested and ready for release, they are merged onto the `develop` branch to await the next release.  The tip of the develop branch always represents the set of features that have been staged for the next release. The version number in all configuration files and the README on the master branch follows [semantic versioning](https://semver.org/) and should always be set to either:

- the current release version, if on `master`. For example, `2.8.5`.
- if on `develop`, the planned next release number with a `beta` designator or release candidate `rc` designator appended as appropriate.  For example, `2.8.6-beta1` or `2.9.0-rc1`.

**feature**. Feature branches should be named as `feature_##_short_title` where `##` is the issue number from github that the feature is implementing.
For example,
`feature_13_r4` might be a new feature being developed independently but intended to be merged into the `develop` branch.

All feature branches should be frequently merged with changes from `develop` to
ensure that the development branch stays up to date with other features that have
been tested and are awaiting release.  Thus, the `develop` branch represents an opportunity
for integration testing of the set of features intended to work together for a
particular release.

## Release process

The release process starts with integration testing in the `develop` branch. Once all
changes that are desired in a release are merged into the `develop` branch, we run
the full set of tests on a clean checkout, and conduct any code reviews needed. A pull request
can be used for that review, or it can be done directly for simpler sets of changes.
Once review is completed, make sure all version numbers are set correctly and 
documentation is complete on `develop`, and then merge to `master`, and tag that merge 
commit as the release. We use tag names of the form `v1.2.3` that reflect the version
of the release.  If we are submitting to CRAN, we generally do not tag the release
until it has been accepted by CRAN.

## Testing

**Unit and integration tests**. We maintain a full suite of unit tests and
integration tests in the `tests` subdirectory.
Any new code developed should include a robust set of unit tests for each public
method, as well as integration tests from new feature sets.  Tests should fully
exercise the feature to ensure that it responds correctly to both good data inputs
as well as various classes of corrupt or bad data.  All tests should pass before
the `develop` branch is merged to master, and all tests should pass before the `master`
branch is tagged as a release.

**Continuous integration**. We can use Travis for some of our repositories.

## Code style

Code should be written to professional standards to enable clean, well-documented,
readable, and maintainable software.  While there has been significant variablility
in the coding styles applied historically, new contributions should strive for
clean code formatting.  Some of the guidelines we follow include:

**Java**. For Java code, follow the [Google Java Style Guide](https://google.github.io/styleguide/javaguide.html), with the single exception that indentation is performed with 4 spaces rather than 2.  When working on a class that
does not follow the conventions, strive to reformat that code module in a single
isolated code commit before starting other code changes.

**Javadoc**. All Java code should be fully documented with JavaDoc comments.  Special
attention should be paid to documentation of the public API for classes.  Documentation
should explain both what the code does, but also why it does it in a particular
way when appropriate.  Class and method documentation should be written to provide
sufficient context for people that are not intimately familiar with the rest of the code.
Class-level documentation often is strengthened through explaining the role of the
class in the architecture.  Avoid using tautological definitions that reuse the name of
a class or method in its definition.  And please be complete.

**R**. Code for R should follow generally accepted R coding conventions.

## Contributor license agreement

In order to clarify the intellectual property license
granted with Contributions from any person or entity, you agree to
a Contributor License Agreement ("CLA") with the Regents of the University of
California (hereafter, the "Regents").

1. Definitions.
   "You" (or "Your") shall mean the copyright owner or legal entity
   authorized by the copyright owner that is making this Agreement
   with the Regents. For legal entities, the entity making a
   Contribution and all other entities that control, are controlled
   by, or are under common control with that entity are considered to
   be a single Contributor. For the purposes of this definition,
   "control" means (i) the power, direct or indirect, to cause the
   direction or management of such entity, whether by contract or
   otherwise, or (ii) ownership of fifty percent (50%) or more of the
   outstanding shares, or (iii) beneficial ownership of such entity.
   "Contribution" shall mean any original work of authorship,
   including any modifications or additions to an existing work, that
   is intentionally submitted by You to the Regents for inclusion
   in, or documentation of, any of the products owned or managed by
   the Regents (the "Work"). For the purposes of this definition,
   "submitted" means any form of electronic, verbal, or written
   communication sent to the Regents or its representatives,
   including but not limited to communication on electronic mailing
   lists, source code control systems, and issue tracking systems that
   are managed by, or on behalf of, the Regents for the purpose of
   discussing and improving the Work, but excluding communication that
   is conspicuously marked or otherwise designated in writing by You
   as "Not a Contribution."
2. Grant of Copyright License. Subject to the terms and conditions of
   this Agreement, You hereby grant to the Regents and to
   recipients of software distributed by the Regents a perpetual,
   worldwide, non-exclusive, no-charge, royalty-free, irrevocable
   copyright license to reproduce, prepare derivative works of,
   publicly display, publicly perform, sublicense, and distribute Your
   Contributions and such derivative works.
3. Grant of Patent License. Subject to the terms and conditions of
   this Agreement, You hereby grant to the Regents and to
   recipients of software distributed by the Regents a perpetual,
   worldwide, non-exclusive, no-charge, royalty-free, irrevocable
   (except as stated in this section) patent license to make, have
   made, use, offer to sell, sell, import, and otherwise transfer the
   Work, where such license applies only to those patent claims
   licensable by You that are necessarily infringed by Your
   Contribution(s) alone or by combination of Your Contribution(s)
   with the Work to which such Contribution(s) was submitted. If any
   entity institutes patent litigation against You or any other entity
   (including a cross-claim or counterclaim in a lawsuit) alleging
   that your Contribution, or the Work to which you have contributed,
   constitutes direct or contributory patent infringement, then any
   patent licenses granted to that entity under this Agreement for
   that Contribution or Work shall terminate as of the date such
   litigation is filed.
4. You represent that you are legally entitled to grant the above
   license. If your employer(s) has rights to intellectual property
   that you create that includes your Contributions, you represent
   that you have received permission to make Contributions on behalf
   of that employer, that your employer has waived such rights for
   your Contributions to the Regents, or that your employer has
   executed a separate Corporate CLA with the Regents.
5. You represent that each of Your Contributions is Your original
   creation (see section 7 for submissions on behalf of others).  You
   represent that Your Contribution submissions include complete
   details of any third-party license or other restriction (including,
   but not limited to, related patents and trademarks) of which you
   are personally aware and which are associated with any part of Your
   Contributions.
