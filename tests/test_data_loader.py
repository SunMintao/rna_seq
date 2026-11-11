def yield_send_demo():
    total = 0
    while True:
        value = yield total  # 接收send来的值
        if value is None:
            break
        total += value
        print(f"累加中间值: {total}")

print("【yield双向通信】")
gen = yield_send_demo()
next(gen)  # 激活生成器

gen.send(10)  # 发送10，打印"累加中间值: 10"
gen.send(20)  # 发送20，打印"累加中间值: 30"
# gen.send(None)  # 结束
