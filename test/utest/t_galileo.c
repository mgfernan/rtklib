#include <assert.h>
#include "../../src/galileo.h"

static void test__compute_sat_pos_from_az_el(void)
{

    struct test_case
    {
        double pos[3];
        double az;
        double el;
        double expected_sat_pos[3];
    } TEST_CASES[] = {
        {
            /* Zenith looking observation */
            .pos = {0.0, 0.0, 0.0},
            .az = 0.0,
            .el = PI / 2,
            .expected_sat_pos = {0.0, 0.0, 25000000.0 - 6378137.0},
        },
        {
            /* North-looking observation */
            .pos = {0.0, 0.0, 0.0},
            .az = 0.0,
            .el = 0.0,
            .expected_sat_pos = {1.312818851305248, 0.0, 25000000.0 - 6378137.0},
        },
        {
            /* South-looking observation */
            .pos = {0.0, 0.0, 0.0},
            .az = PI,
            .el = 0.0,
            .expected_sat_pos = {-1.312818851305248, 0.0, 25000000.0 - 6378137.0},
        },
        {
            /* East-looking observation */
            .pos = {0.0, 0.0, 0.0},
            .az = PI/2.0,
            .el = 0.0,
            .expected_sat_pos = {0.0, 1.312818851305248, 25000000.0 - 6378137.0},
        },
        {
            /* West-looking observation */
            .pos = {0.0, 0.0, 0.0},
            .az = 3.0 * PI / 2.0,
            .el = 0.0,
            .expected_sat_pos = {0.0, -1.312818851305248, 25000000.0 - 6378137.0},
        },
    };

    int n_test_cases = sizeof(TEST_CASES) / sizeof(TEST_CASES[0]);

    for (int i = 0; i < n_test_cases; i ++) {
        struct test_case *tc = &TEST_CASES[i];
        double sat_pos[3];
        compute_sat_pos_from_az_el(tc->pos, tc->az, tc->el, sat_pos);

        assert(sat_pos[0] - tc->expected_sat_pos[0] < 1e-6);
        assert(sat_pos[1] - tc->expected_sat_pos[1] < 1e-6);
        assert(sat_pos[2] - tc->expected_sat_pos[2] < 1e-6);
    }

}


int main(void)
{
    test__compute_sat_pos_from_az_el();

    return 0;
}