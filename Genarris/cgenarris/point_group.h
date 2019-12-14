typedef struct
{
	char point_group[6];
	char sub_groups[32][6];
	int len_sub_groups;
}point_group_database;

static const point_group_database point_group_dbs[] = {
{
        .point_group = "m-3m\0",
         .len_sub_groups = 25 ,
        .sub_groups = { "m-3m\0", "-43m\0", "432\0", "m-3\0", "4/mmm\0", "23\0", "-3m\0", "-42m\0", "4mm\0", "422\0", "4/m\0", "mmm\0", "3m\0", "32\0", "-3\0", "-4\0", "mm2\0", "222\0", "4\0", "2/m\0", "3\0", "m\0", "2\0", "-1\0", "1\0"}
 },
{
        .point_group = "-43m\0",
         .len_sub_groups = 11 ,
        .sub_groups = { "-43m\0", "23\0", "-42m\0", "3m\0", "-4\0", "mm2\0", "222\0", "3\0", "m\0", "2\0", "1\0"}
 },
{
        .point_group = "432\0",
         .len_sub_groups = 9 ,
        .sub_groups = { "432\0", "23\0", "422\0", "32\0", "222\0", "4\0", "3\0", "2\0", "1\0"}
 },
{
        .point_group = "m-3\0",
         .len_sub_groups = 12 ,
        .sub_groups = { "m-3\0", "23\0", "mmm\0", "-3\0", "mm2\0", "222\0", "2/m\0", "3\0", "m\0", "2\0", "-1\0", "1\0"}
 },
{
        .point_group = "6/mmm\0",
         .len_sub_groups = 20 ,
        .sub_groups = { "6/mmm\0", "6mm\0", "-6m2\0", "622\0", "6/m\0", "-3m\0", "mmm\0", "3m\0", "32\0", "-3\0", "-6\0", "6\0", "mm2\0", "222\0", "2/m\0", "3\0", "m\0", "2\0", "-1\0", "1\0"}
 },
{
        .point_group = "4/mmm\0",
         .len_sub_groups = 15 ,
        .sub_groups = { "4/mmm\0", "-42m\0", "4mm\0", "422\0", "4/m\0", "mmm\0", "-4\0", "mm2\0", "222\0", "4\0", "2/m\0", "m\0", "2\0", "-1\0", "1\0"}
 },
{
        .point_group = "23\0",
         .len_sub_groups = 5 ,
        .sub_groups = { "23\0", "222\0", "3\0", "2\0", "1\0"}
 },
{
        .point_group = "6mm\0",
         .len_sub_groups = 9 ,
        .sub_groups = { "6mm\0", "-3\0", "6\0", "2/m\0", "3\0", "m\0", "2\0", "-1\0", "1\0"}
 },
{
        .point_group = "-6m2\0",
         .len_sub_groups = 9 ,
        .sub_groups = { "-6m2\0", "3m\0", "32\0", "-6\0", "mm2\0", "3\0", "m\0", "2\0", "1\0"}
 },
{
        .point_group = "622\0",
         .len_sub_groups = 7 ,
        .sub_groups = { "622\0", "32\0", "6\0", "222\0", "3\0", "2\0", "1\0"}
 },
{
        .point_group = "6/m\0",
         .len_sub_groups = 10 ,
        .sub_groups = { "6/m\0", "-3\0", "-6\0", "6\0", "2/m\0", "3\0", "m\0", "2\0", "-1\0", "1\0"}
 },
{
        .point_group = "-3m\0",
         .len_sub_groups = 10 ,
        .sub_groups = { "-3m\0", "3m\0", "32\0", "-3\0", "2/m\0", "3\0", "m\0", "2\0", "-1\0", "1\0"}
 },
{
        .point_group = "-42m\0",
         .len_sub_groups = 7 ,
        .sub_groups = { "-42m\0", "-4\0", "mm2\0", "222\0", "m\0", "2\0", "1\0"}
 },
{
        .point_group = "4mm\0",
         .len_sub_groups = 6 ,
        .sub_groups = { "4mm\0", "mm2\0", "4\0", "m\0", "2\0", "1\0"}
 },
{
        .point_group = "422\0",
         .len_sub_groups = 5 ,
        .sub_groups = { "422\0", "222\0", "4\0", "2\0", "1\0"}
 },
{
        .point_group = "4/m\0",
         .len_sub_groups = 8 ,
        .sub_groups = { "4/m\0", "-4\0", "4\0", "2/m\0", "m\0", "2\0", "-1\0", "1\0"}
 },
{
        .point_group = "mmm\0",
         .len_sub_groups = 8 ,
        .sub_groups = { "mmm\0", "mm2\0", "222\0", "2/m\0", "m\0", "2\0", "-1\0", "1\0"}
 },
{
        .point_group = "3m\0",
         .len_sub_groups = 4 ,
        .sub_groups = { "3m\0", "3\0", "m\0", "1\0"}
 },
{
        .point_group = "32\0",
         .len_sub_groups = 4 ,
        .sub_groups = { "32\0", "3\0", "2\0", "1\0"}
 },
{
        .point_group = "-3\0",
         .len_sub_groups = 4 ,
        .sub_groups = { "-3\0", "3\0", "-1\0", "1\0"}
 },
{
        .point_group = "-6\0",
         .len_sub_groups = 4 ,
        .sub_groups = { "-6\0", "3\0", "m\0", "1\0"}
 },
{
        .point_group = "6\0",
         .len_sub_groups = 4 ,
        .sub_groups = { "6\0", "3\0", "2\0", "1\0"}
 },
{
        .point_group = "-4\0",
         .len_sub_groups = 3 ,
        .sub_groups = { "-4\0", "2\0", "1\0"}
 },
{
        .point_group = "mm2\0",
         .len_sub_groups = 4 ,
        .sub_groups = { "mm2\0", "m\0", "2\0", "1\0"}
 },
{
        .point_group = "222\0",
         .len_sub_groups = 3 ,
        .sub_groups = { "222\0", "2\0", "1\0"}
 },
{
        .point_group = "4\0",
         .len_sub_groups = 3 ,
        .sub_groups = { "4\0", "2\0", "1\0"}
 },
{
        .point_group = "2/m\0",
         .len_sub_groups = 5 ,
        .sub_groups = { "2/m\0", "m\0", "2\0", "-1\0", "1\0"}
 },
{
        .point_group = "3\0",
         .len_sub_groups = 2 ,
        .sub_groups = { "3\0", "1\0"}
 },
{
        .point_group = "m\0",
         .len_sub_groups = 2 ,
        .sub_groups = { "m\0", "1\0"}
 },
{
        .point_group = "2\0",
         .len_sub_groups = 2 ,
        .sub_groups = { "2\0", "1\0"}
 },
{
        .point_group = "-1\0",
         .len_sub_groups = 2 ,
        .sub_groups = { "-1\0", "1\0"}
 },
{
        .point_group = "1\0",
         .len_sub_groups = 1 ,
        .sub_groups = { "1\0"}
 }
};
