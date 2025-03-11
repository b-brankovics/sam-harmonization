use 5.026001;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'biointsam',
    VERSION_FROM      => 'lib/biointsam.pm', # finds $VERSION, requires EU::MM from perl >= 5.5
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ABSTRACT_FROM     => 'lib/biointsam.pm', # retrieve abstract from module
    AUTHOR            => 'b-brankovics <6728856+b-brankovics@users.noreply.github.com>',
#    EXE_FILES         => [ "sam-score.pl" ], #[ map "bin/$_", qw(sam-ref-plot.pl) ],
    EXE_FILES         => [ glob("*.pl") ],
    LICENSE           => 'gpl',
    #Value must be from legacy list of licenses here
    #http://search.cpan.org/perldoc?Module%3A%3ABuild%3A%3AAPI
);
