function test_noise_seed()
%TEST_NOISE_SEED nmpc_base owns the measurement-noise seed.
%   The scripts used to call rng(seed) and then nmpc_base, which drew with
%   randn from the global stream. nmpc_base now draws from a stream of its own;
%   this test checks that the realisation is the one rng(seed) gave, that it
%   does not depend on what ran before, and that the global stream is left
%   alone.

    repo = fileparts(fileparts(fileparts(mfilename("fullpath"))));
    addpath(genpath(repo));

    calls = {
        "campaign",  {}
        "noiseless", {"sigma_y", [0 0 0]}
        "schedule",  {"tf", 30, "set_setpoint", false, "load_lqr", false}
        "seed 345",  {"noise_seed", 345}
    };
    for k = 1:size(calls, 1)
        args = calls{k, 2};
        b = nmpc_base(args{:});
        rng(b.noise_seed);
        expected = randn(b.N, b.nx) .* b.sigma_y;
        assert(isequal(b.noise, expected), "%s: noise is not the rng(%d) realisation", ...
            calls{k, 1}, b.noise_seed);
        fprintf("%s: PASS (seed %d)\n", calls{k, 1}, b.noise_seed);
    end

    % Two bases in one session get the same realisation, whatever ran between.
    rng(999);
    first = nmpc_base(sigma_y = [0 0 0]); %#ok<NASGU>
    randn(50, 1);
    second = nmpc_base();
    rng(123);
    assert(isequal(second.noise, randn(second.N, second.nx) .* second.sigma_y), ...
        "a second base in one session got a different realisation");
    fprintf("second base in a session: PASS\n");

    rng(999);
    nmpc_base();
    assert(rng().Seed == 999, "nmpc_base changed the global random state");
    fprintf("global stream untouched: PASS\n");

    % No script seeds the noise itself any more.
    scripts = dir(fullfile(repo, "*.m"));
    seeded = string({scripts.name});
    seeded = seeded(arrayfun(@(f) contains(fileread(fullfile(repo, f)), "rng(123)"), seeded));
    assert(isempty(seeded), "these scripts still call rng(123): %s", strjoin(seeded, ", "));
    fprintf("no rng(123) in the entry scripts: PASS\n");
end
