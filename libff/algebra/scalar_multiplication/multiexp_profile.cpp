#include <cstdio>
#include <vector>

#include <libff/algebra/curves/bn128/bn128_pp.hpp>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/rng.hpp>

#define INSTANCE_SIZE 10000
#define INSTANCE_COUNT 10

using namespace libff;

template<typename GroupT, typename FieldT>
void profile_multiexp()
{
    std::vector<std::vector<GroupT> > group_elements(INSTANCE_COUNT);
    std::vector<std::vector<FieldT> > scalars(INSTANCE_COUNT);

    enter_block("Generating test data");
    for (int i = 0; i < INSTANCE_COUNT; i++) {
        GroupT x = GroupT::random_element();
        for (int j = 0; j < INSTANCE_SIZE; j++) {
            group_elements[i].push_back(x);
            scalars[i].push_back(SHA512_rng<FieldT>(i * INSTANCE_SIZE + j));
        }
    }
    leave_block("Generating test data");


    enter_block("Running naive_exp");
    std::vector<GroupT> naive_answers;
    for (int i = 0; i < INSTANCE_COUNT; i++) {
        naive_answers.push_back(multi_exp<GroupT, FieldT>(
            group_elements[i].cbegin(), group_elements[i].cend(),
            scalars[i].cbegin(), scalars[i].cend(),
            1, false));
    }
    leave_block("Running naive_exp");

    enter_block("Running multi_exp");
    std::vector<GroupT> multiexp_answers;
    for (int i = 0; i < INSTANCE_COUNT; i++) {
        multiexp_answers.push_back(multi_exp<GroupT, FieldT>(
            group_elements[i].cbegin(), group_elements[i].cend(),
            scalars[i].cbegin(), scalars[i].cend(),
            1, true));
    }
    leave_block("Running multi_exp");

    bool answers_correct = true;
    for (int i = 0; i < INSTANCE_COUNT; i++) {
        answers_correct = answers_correct &&
            (naive_answers[i] == multiexp_answers[i]);
    }

    if (answers_correct) {
        printf("Answers correct\n");
    } else {
        printf("!!!!!!!!!!!!!!!!!\n");
        printf("Answers INCORRECT\n");
        printf("!!!!!!!!!!!!!!!!!\n");
    }
}

int main(void)
{
    print_compilation_info();
    start_profiling();

    printf("Using %d instances of size %d each\n", INSTANCE_COUNT, INSTANCE_SIZE);

    printf("Profiling BN128_G1\n");
    bn128_pp::init_public_params();
    profile_multiexp<G1<bn128_pp>, Fr<bn128_pp> >();

    return 0;
}
