import itertools

all_ops_cat = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9"]
cross_terms_ = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_ops_cat, 2)]
print("Defined cross terms:", cross_terms_)

Cross_terms_FM= ['FM0vsFM1', 'FM0vsFM2', 'FM0vsFM3', 'FM0vsFM4', 'FM0vsFM5', 'FM0vsFM7', 'FM0vsFM8', 'FM0vsFM9'
                 ,'FM1vsFM2', 'FM1vsFM3', 'FM1vsFM4', 'FM1vsFM5', 'FM1vsFM7', 'FM1vsFM8', 'FM1vsFM9',
                 'FM2vsFM3', 'FM2vsFM4', 'FM2vsFM5', 'FM2vsFM7', 'FM2vsFM8', 'FM2vsFM9',
                 'FM3vsFM4', 'FM3vsFM5', 'FM3vsFM7', 'FM3vsFM8', 'FM3vsFM9',
                 'FM4vsFM5', 'FM4vsFM7', 'FM4vsFM8', 'FM4vsFM9',
                 'FM5vsFM7', 'FM5vsFM8', 'FM5vsFM9', 'FM7vsFM8', 'FM7vsFM9', 'FM8vsFM9']

cross_term_FM=  ['FM0vsFM1', 'FM0vsFM2', 'FM0vsFM3', 'FM0vsFM4', 'FM0vsFM5', 'FM0vsFM7', 'FM0vsFM8', 'FM0vsFM9', 'FM1vsFM2', 'FM1vsFM3', 'FM1vsFM4', 'FM1vsFM5', 'FM1vsFM7', 'FM1vsFM8', 'FM1vsFM9', 'FM2vsFM3', 'FM2vsFM4', 'FM2vsFM5', 'FM2vsFM7', 'FM2vsFM8', 'FM2vsFM9', 'FM3vsFM4', 'FM3vsFM5', 'FM3vsFM7', 'FM3vsFM8', 'FM3vsFM9', 'FM4vsFM5', 'FM4vsFM7', 'FM4vsFM8', 'FM4vsFM9', 'FM5vsFM7', 'FM5vsFM8', 'FM5vsFM9', 'FM7vsFM8', 'FM7vsFM9', 'FM8vsFM9']

Cross_terms_FS= ['FS0vsFS1', 'FS0vsFS2', 'FS1vsFS2']

Cross_terms_FT= ['FT0vsFT1', 'FT0vsFT2', 'FT0vsFT3', 'FT0vsFT4', 'FT0vsFT5', 'FT0vsFT6',
                    'FT1vsFT2', 'FT1vsFT3', 'FT1vsFT4', 'FT1vsFT5', 'FT1vsFT6',
                    'FT2vsFT3', 'FT2vsFT4', 'FT2vsFT5', 'FT2vsFT6',
                    'FT3vsFT4', 'FT3vsFT5', 'FT3vsFT6',
                    'FT4vsFT5', 'FT4vsFT6',
                    'FT5vsFT6']