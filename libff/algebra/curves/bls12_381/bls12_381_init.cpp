#include <libff/algebra/curves/bls12_381/bls12_381_g1.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_g2.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_init.hpp>

namespace libff {

bls12_381_Fq bls12_381_coeff_b;
bls12_381_Fq2 bls12_381_twist;
bls12_381_Fq2 bls12_381_twist_coeff_b;
bls12_381_Fq bls12_381_twist_mul_by_b_c0;
bls12_381_Fq bls12_381_twist_mul_by_b_c1;
bls12_381_Fq2 bls12_381_twist_mul_by_q_X;
bls12_381_Fq2 bls12_381_twist_mul_by_q_Y;

bigint<bls12_381_q_limbs> bls12_381_ate_loop_count;
bool bls12_381_ate_is_loop_count_neg;
bigint<12*bls12_381_q_limbs> bls12_381_final_exponent;
bigint<bls12_381_q_limbs> bls12_381_final_exponent_z;
bool bls12_381_final_exponent_is_z_neg;

void init_bls12_381_params()
{
    init_bls12_381_fields();

    /* choice of short Weierstrass curve and its twist */

    bls12_381_coeff_b = bls12_381_Fq("4");
    bls12_381_twist = bls12_381_Fq2(bls12_381_Fq("1"), bls12_381_Fq("1"));
    bls12_381_twist_coeff_b = bls12_381_coeff_b * bls12_381_twist;
    bls12_381_twist_mul_by_b_c0 = bls12_381_coeff_b * bls12_381_Fq2::non_residue;
    bls12_381_twist_mul_by_b_c1 = bls12_381_coeff_b * bls12_381_Fq2::non_residue;
    bls12_381_twist_mul_by_q_X = bls12_381_Fq2(bls12_381_Fq("0"),
                                               bls12_381_Fq("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939437"));
    bls12_381_twist_mul_by_q_Y = bls12_381_Fq2(bls12_381_Fq("2973677408986561043442465346520108879172042883009249989176415018091420807192182638567116318576472649347015917690530"),
                                               bls12_381_Fq("1028732146235106349975324479215795277384839936929757896155643118032610843298655225875571310552543014690878354869257"));


    /* choice of group G1 */
    bls12_381_G1::G1_zero = bls12_381_G1(bls12_381_Fq::zero(),
                                     bls12_381_Fq::one(),
                                     bls12_381_Fq::zero());
    bls12_381_G1::G1_one = bls12_381_G1(bls12_381_Fq("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507"),
                                    bls12_381_Fq("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569"),
                                    bls12_381_Fq::one());
    // Cofactor
    bls12_381_G1::h = bigint<bls12_381_G1::h_limbs>("76329603384216526031706109802092473003");


    // TODO: wNAF window table
    bls12_381_G1::wnaf_window_table.resize(0);
    bls12_381_G1::wnaf_window_table.push_back(11);
    bls12_381_G1::wnaf_window_table.push_back(24);
    bls12_381_G1::wnaf_window_table.push_back(60);
    bls12_381_G1::wnaf_window_table.push_back(127);

    // TODO: fixed-base exponentiation table
    bls12_381_G1::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 4.99]
    bls12_381_G1::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [4.99, 10.99]
    bls12_381_G1::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.99, 32.29]
    bls12_381_G1::fixed_base_exp_window_table.push_back(11);
    // window 4 is unbeaten in [32.29, 55.23]
    bls12_381_G1::fixed_base_exp_window_table.push_back(32);
    // window 5 is unbeaten in [55.23, 162.03]
    bls12_381_G1::fixed_base_exp_window_table.push_back(55);
    // window 6 is unbeaten in [162.03, 360.15]
    bls12_381_G1::fixed_base_exp_window_table.push_back(162);
    // window 7 is unbeaten in [360.15, 815.44]
    bls12_381_G1::fixed_base_exp_window_table.push_back(360);
    // window 8 is unbeaten in [815.44, 2373.07]
    bls12_381_G1::fixed_base_exp_window_table.push_back(815);
    // window 9 is unbeaten in [2373.07, 6977.75]
    bls12_381_G1::fixed_base_exp_window_table.push_back(2373);
    // window 10 is unbeaten in [6977.75, 7122.23]
    bls12_381_G1::fixed_base_exp_window_table.push_back(6978);
    // window 11 is unbeaten in [7122.23, 57818.46]
    bls12_381_G1::fixed_base_exp_window_table.push_back(7122);
    // window 12 is never the best
    bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // window 13 is unbeaten in [57818.46, 169679.14]
    bls12_381_G1::fixed_base_exp_window_table.push_back(57818);
    // window 14 is never the best
    bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // window 15 is unbeaten in [169679.14, 439758.91]
    bls12_381_G1::fixed_base_exp_window_table.push_back(169679);
    // window 16 is unbeaten in [439758.91, 936073.41]
    bls12_381_G1::fixed_base_exp_window_table.push_back(439759);
    // window 17 is unbeaten in [936073.41, 4666554.74]
    bls12_381_G1::fixed_base_exp_window_table.push_back(936073);
    // window 18 is never the best
    bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4666554.74, 7580404.42]
    bls12_381_G1::fixed_base_exp_window_table.push_back(4666555);
    // window 20 is unbeaten in [7580404.42, 34552892.20]
    bls12_381_G1::fixed_base_exp_window_table.push_back(7580404);
    // window 21 is never the best
    bls12_381_G1::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [34552892.20, inf]
    bls12_381_G1::fixed_base_exp_window_table.push_back(34552892);


    /* choice of group G2 */
    bls12_381_G2::G2_zero = bls12_381_G2(bls12_381_Fq2::zero(),
                                         bls12_381_Fq2::one(),
                                         bls12_381_Fq2::zero());

    // simple G2 generator
    bls12_381_G2::G2_one = bls12_381_G2(bls12_381_Fq2(bls12_381_Fq("352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160"),
                                                      bls12_381_Fq("3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758")),
                                        bls12_381_Fq2(bls12_381_Fq("1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905"),
                                                      bls12_381_Fq("927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582")),
                                        bls12_381_Fq2::one());
    // Cofactor
    bls12_381_G2::h = bigint<bls12_381_G2::h_limbs>("305502333931268344200999753193121504214466019254188142667664032982267604182971884026507427359259977847832272839041616661285803823378372096355777062779109");


    // TODO: wNAF window table
    bls12_381_G2::wnaf_window_table.resize(0);
    bls12_381_G2::wnaf_window_table.push_back(5);
    bls12_381_G2::wnaf_window_table.push_back(15);
    bls12_381_G2::wnaf_window_table.push_back(39);
    bls12_381_G2::wnaf_window_table.push_back(109);

    // TODO: fixed-base exponentiation table
    bls12_381_G2::fixed_base_exp_window_table.resize(0);
    // window 1 is unbeaten in [-inf, 5.10]
    bls12_381_G2::fixed_base_exp_window_table.push_back(1);
    // window 2 is unbeaten in [5.10, 10.43]
    bls12_381_G2::fixed_base_exp_window_table.push_back(5);
    // window 3 is unbeaten in [10.43, 25.28]
    bls12_381_G2::fixed_base_exp_window_table.push_back(10);
    // window 4 is unbeaten in [25.28, 59.00]
    bls12_381_G2::fixed_base_exp_window_table.push_back(25);
    // window 5 is unbeaten in [59.00, 154.03]
    bls12_381_G2::fixed_base_exp_window_table.push_back(59);
    // window 6 is unbeaten in [154.03, 334.25]
    bls12_381_G2::fixed_base_exp_window_table.push_back(154);
    // window 7 is unbeaten in [334.25, 742.58]
    bls12_381_G2::fixed_base_exp_window_table.push_back(334);
    // window 8 is unbeaten in [742.58, 2034.40]
    bls12_381_G2::fixed_base_exp_window_table.push_back(743);
    // window 9 is unbeaten in [2034.40, 4987.56]
    bls12_381_G2::fixed_base_exp_window_table.push_back(2034);
    // window 10 is unbeaten in [4987.56, 8888.27]
    bls12_381_G2::fixed_base_exp_window_table.push_back(4988);
    // window 11 is unbeaten in [8888.27, 26271.13]
    bls12_381_G2::fixed_base_exp_window_table.push_back(8888);
    // window 12 is unbeaten in [26271.13, 39768.20]
    bls12_381_G2::fixed_base_exp_window_table.push_back(26271);
    // window 13 is unbeaten in [39768.20, 106275.75]
    bls12_381_G2::fixed_base_exp_window_table.push_back(39768);
    // window 14 is unbeaten in [106275.75, 141703.40]
    bls12_381_G2::fixed_base_exp_window_table.push_back(106276);
    // window 15 is unbeaten in [141703.40, 462422.97]
    bls12_381_G2::fixed_base_exp_window_table.push_back(141703);
    // window 16 is unbeaten in [462422.97, 926871.84]
    bls12_381_G2::fixed_base_exp_window_table.push_back(462423);
    // window 17 is unbeaten in [926871.84, 4873049.17]
    bls12_381_G2::fixed_base_exp_window_table.push_back(926872);
    // window 18 is never the best
    bls12_381_G2::fixed_base_exp_window_table.push_back(0);
    // window 19 is unbeaten in [4873049.17, 5706707.88]
    bls12_381_G2::fixed_base_exp_window_table.push_back(4873049);
    // window 20 is unbeaten in [5706707.88, 31673814.95]
    bls12_381_G2::fixed_base_exp_window_table.push_back(5706708);
    // window 21 is never the best
    bls12_381_G2::fixed_base_exp_window_table.push_back(0);
    // window 22 is unbeaten in [31673814.95, inf]
    bls12_381_G2::fixed_base_exp_window_table.push_back(31673815);



    /* pairing parameters */

    bls12_381_ate_loop_count = bigint<bls12_381_q_limbs>("15132376222941642752");
    bls12_381_ate_is_loop_count_neg = true;
    bls12_381_final_exponent = bigint<12*bls12_381_q_limbs>("322277361516934140462891564586510139908379969514828494218366688025288661041104682794998680497580008899973249814104447692778988208376779573819485263026159588510513834876303014016798809919343532899164848730280942609956670917565618115867287399623286813270357901731510188149934363360381614501334086825442271920079363289954510565375378443704372994881406797882676971082200626541916413184642520269678897559532260949334760604962086348898118982248842634379637598665468817769075878555493752214492790122785850202957575200176084204422751485957336465472324810982833638490904279282696134323072515220044451592646885410572234451732790590013479358343841220074174848221722017083597872017638514103174122784843925578370430843522959600095676285723737049438346544753168912974976791528535276317256904336520179281145394686565050419250614107803233314658825463117900250701199181529205942363159325765991819433914303908860460720581408201373164047773794825411011922305820065611121544561808414055302212057471395719432072209245600258134364584636810093520285711072578721435517884103526483832733289802426157301542744476740008494780363354305116978805620671467071400711358839553375340724899735460480144599782014906586543813292157922220645089192130209334926661588737007768565838519456601560804957985667880395221049249803753582637708560");
    bls12_381_final_exponent_z = bigint<bls12_381_q_limbs>("15132376222941642752");
    bls12_381_final_exponent_is_z_neg = true;
}
} // namespace libff
